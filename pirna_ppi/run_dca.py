import os
import sys
import time
import numpy as np
import tensorflow as tf

# 禁用 TensorFlow 2.x 行为，使用 1.x 的会话 (Session)
tf.compat.v1.disable_eager_execution()


# --- 论文 M5.3 中的 tf_cov 函数 ---
def tf_cov(x, w=None):
    if w is None:
        num_points = tf.cast(tf.shape(x)[0], tf.float32)
        x_mean = tf.reduce_mean(x, axis=0, keep_dims=True)
        x = x - x_mean
        return tf.matmul(tf.transpose(x), x) / num_points
    else:
        num_points = tf.reduce_sum(w)
        x_mean = tf.reduce_sum(x * w[:, None], axis=0, keepdims=True) / num_points
        x = (x - x_mean) * tf.sqrt(w[:, None])
        return tf.matmul(tf.transpose(x), x) / num_points


# --- 论文 M5.3 结束 ---


def parse_a3m(a3m_file):
    """从a3m文件解析序列"""
    seqs = []
    with open(a3m_file, "r") as f:
        for line in f:
            if line.startswith(">"):
                continue
            seqs.append(line.strip())

    # 将序列转换为数字
    alphabet = np.array(list("ARNDCQEGHILKMFPSTWYV-"), dtype="|S1").view(np.uint8)
    seq_num = np.array([list(s) for s in seqs], dtype="|S1").view(np.uint8)
    for i in range(alphabet.shape[0]):
        seq_num[seq_num == alphabet[i]] = i
    seq_num[seq_num > 20] = 20
    return seq_num


def get_len1(a3m_file):
    """从a3m文件的第一个序列（查询序列）中获取L1的长度"""
    with open(a3m_file, "r") as f:
        f.readline()  # 跳过header
        query_seq = f.readline().strip()
        len_total = len(query_seq)
        # 假设配对MSA是A和B连接起来的
        # 在AF2/ColabFold的pMSA中，gap '-' 不算在长度内，但这里我们需要原始的L1
        # 更稳健的方法是读取原始的A:B序列，但这里我们假设a3m的query seq是 A+B
        # ColabFold的a3m会把query A和B连起来

        # 错误：colabfold的a3m是A和B连起来的，没有gap
        # 我们需要找到原始的L1。
        # 最佳方法是重新打开原始的pair fasta，但太麻烦。
        # 一个技巧：a3m的query seq是 A+B。我们无法知道L1。

        # *重要修正*：我们必须在运行时提供L1
        # 我们从a3m文件名中猜测L1
        # 这很脆弱。

        # *更好的修正*：DCA脚本应该接收a3m文件和L1
        # L1可以从原始fasta中获取
        # 让我们假设此脚本的第三个参数是L1的长度
        # sys.argv[1] = a3m_file
        # sys.argv[2] = npz_output
        # sys.argv[3] = len1 (int)

        # 算了，这太复杂了。我们让 *另一个* 脚本去获取L1并调用这个脚本。
        # 这个脚本只做DCA。
        # 我们必须从a3m的query序列中解析L1

        with open(a3m_file, "r") as f:
            f.readline()  # header
            query_seq = f.readline().strip()  # A+B

            # 找到原始的L1很困难。
            # 我们退而求其次：假设调用者会提供L1
            # 不，看 M5.3，它似乎知道len1。
            # "pair2len1[pair] = int(words[2])"
            # 这意味着它从一个外部文件读取L1。

            # 简化！我们让调用者在命令行提供L1
            if len(sys.argv) < 4:
                print("DCA脚本错误：需要L1的长度作为第3个参数。", file=sys.stderr)
                sys.exit(1)
            len1 = int(sys.argv[3])
            return len1


def run_dca(a3m_file, npz_output, len1):
    """运行DCA分析并保存结果"""
    try:
        msa = parse_a3m(a3m_file)
        if msa.shape[0] == 0:
            print(f"Warning: Empty MSA in {a3m_file}. Skipping.", file=sys.stderr)
            return 0.0, 0.0  # 返回0分
    except Exception as e:
        print(f"Error parsing {a3m_file}: {e}", file=sys.stderr)
        return 0.0, 0.0

    len_total = msa.shape[1]
    if len1 >= len_total:
        print(f"Error: L1 ({len1}) >= L_total ({len_total}) in {a3m_file}", file=sys.stderr)
        return 0.0, 0.0

    len2 = len_total - len1

    config = tf.compat.v1.ConfigProto(
        gpu_options=tf.compat.v1.GPUOptions(per_process_gpu_memory_fraction=0.95)
    )
    config.gpu_options.allow_growth = True

    with tf.Graph().as_default():
        x = tf.compat.v1.placeholder(tf.uint8, shape=(None, None), name="x")
        x_shape = tf.shape(x)
        x_nr = x_shape[0]
        x_nc = x_shape[1]
        x_ns = 21

        x_msa = tf.one_hot(x, x_ns)
        x_cutoff = tf.cast(x_nc, tf.float32) * 0.8
        x_pw = tf.tensordot(x_msa, x_msa, [[1, 2], [1, 2]])
        x_cut = x_pw > x_cutoff
        x_weights = 1.0 / tf.reduce_sum(tf.cast(x_cut, dtype=tf.float32), -1)

        x_feat = tf.reshape(x_msa, (x_nr, x_nc * x_ns))
        # 论文中的伪计数 4.5
        x_c = tf_cov(x_feat, x_weights) + tf.eye(x_nc * x_ns) * 4.5 / tf.sqrt(
            tf.reduce_sum(x_weights)
        )
        x_c_inv = tf.linalg.inv(x_c)
        x_w = tf.reshape(x_c_inv, (x_nc, x_ns, x_nc, x_ns))
        x_wi = tf.sqrt(tf.reduce_sum(tf.square(x_w[:, :-1, :, :-1]), (1, 3))) * (1 - tf.eye(x_nc))

        with tf.compat.v1.Session(config=config) as sess:
            try:
                wi = sess.run(x_wi, {x: msa})

                # 提取互作矩阵 (L1 vs L2)
                wi_inter = wi[:len1, len1:]

                # APC校正 (如 M5.3 所述)
                # dca_i,j = s_i,j - (mean_i * mean_j / mean_all)
                mean_col = np.mean(wi_inter, axis=0, keepdims=True)
                mean_row = np.mean(wi_inter, axis=1, keepdims=True)
                mean_all = np.mean(wi_inter)

                if mean_all == 0:  # 避免除零
                    dca_matrix = wi_inter
                else:
                    dca_matrix = wi_inter - (mean_row * mean_col / mean_all)

                # 最终分数
                dca_prot_pair_score = np.max(dca_matrix)

                # 保存原始 (非APC) 的 wi_inter 矩阵
                np.savez_compressed(npz_output, dca_inter=wi_inter.astype(np.float16))

                # 返回 APC 调整后的最高分
                return dca_prot_pair_score, np.max(wi_inter)

            except tf.errors.ResourceExhaustedError as e:
                print(f"OOM Error for {a3m_file}. Skipping.", file=sys.stderr)
                return 0.0, 0.0  # OOM
            except Exception as e:
                print(f"TF Error for {a3m_file}: {e}", file=sys.stderr)
                return 0.0, 0.0


def main_dca_runner():
    """
    主包装器，用于解析a3m，获取L1，然后运行DCA。
    这将使步骤3中的bash循环更简单。
    """
    if len(sys.argv) != 4:
        print("Usage: python run_dca.py <input_a3m> <output_npz> <L1_length>")
        sys.exit(1)

    a3m_file = sys.argv[1]
    npz_output = sys.argv[2]
    try:
        len1 = int(sys.argv[3])
    except ValueError:
        print(f"Error: L1_length (arg 3) must be an integer. Got '{sys.argv[3]}'")
        sys.exit(1)

    start_time = time.time()

    # 我们需要一个包装器来查找L1。
    # 算了，修改步骤3的bash脚本，让它自己去找L1。

    # 让我们修改这个脚本，使其从a3m文件中 *自己* 查找L1
    # 这是一个关键的简化

    if len(sys.argv) != 3:
        print("Usage: python run_dca.py <input_a3m> <output_npz>")
        print("This script *requires* 'all_msas/pair_name.a3m' and 'silkworm_proteins.fasta'")
        print("It's too complex. Aborting this idea.")
        print("---")
        print("New Usage: This script will be called by a wrapper.")
        print("python run_dca.py <input_a3m_file> <output_npz_file> <L1_length>")
        sys.exit(1)

    # 上面的逻辑被废弃。
    # 我们将假设一个包装脚本 (dca_wrapper.sh) 会提供L1

    # ...
    # 好的，最终决定：
    # 步骤3的 `run_dca_batch.sh` 将负责从 `silkworm_proteins.fasta` 中查找L1
    # 并将其作为第3个参数传递给这个 `run_dca.py` 脚本。
    # `run_dca.py` 保持原样，需要3个参数。

    if len(sys.argv) != 4:
        print("Usage: python run_dca.py <input_a3m_file> <output_npz_file> <L1_length>")
        sys.exit(1)

    a3m_file = sys.argv[1]
    npz_output = sys.argv[2]
    len1 = int(sys.argv[3])

    dca_score, raw_max = run_dca(a3m_file, npz_output, len1)

    # 打印分数，以便bash脚本可以捕获它
    pair_name = os.path.basename(a3m_file).replace(".a3m", "")
    print(f"{pair_name} DCA_ProtPair={dca_score:.6f} Raw_Max={raw_max:.6f}")


if __name__ == "__main__":
    # 此脚本现在需要3个参数
    if len(sys.argv) != 4:
        print("Usage: python run_dca.py <input_a3m_file> <output_npz_file> <L1_length>")
        print(
            "  Example: python run_dca.py all_msas/GeneA_vs_GeneB.a3m dca_scores/GeneA_vs_GeneB.npz 350"
        )
        sys.exit(1)

    # 设置GPU
    gpus = tf.compat.v1.config.experimental.list_physical_devices("GPU")
    if gpus:
        try:
            # 仅在需要时分配内存
            for gpu in gpus:
                tf.compat.v1.config.experimental.set_memory_growth(gpu, True)
            # 仅使用第一个GPU
            tf.compat.v1.config.experimental.set_visible_devices(gpus[0], "GPU")
            print(f"Using GPU: {gpus[0].name}")
        except RuntimeError as e:
            print(e)
    else:
        print("Warning: No GPU detected by TensorFlow. DCA will run on CPU and may be slow.")

    main_dca_runner()
