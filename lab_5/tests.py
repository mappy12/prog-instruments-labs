import math
import logging

from scipy.special import gammainc

from consts import *

logging.basicConfig(
    level=logging.DEBUG,
    format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
    handlers=[
        logging.FileHandler("tests.log", encoding="utf-8"),
        logging.StreamHandler()
    ]
)

logger = logging.getLogger(__name__)

def read_file(filename: str) -> str:
    """
    Reads the sequence

    :param filename: Path to the file to read.
    :return: The sequence
    """

    logger.info(f"Reading file: {filename}")

    try:

        with open(filename, 'r', encoding='utf-8') as file:
            data = file.read()
            logger.debug(
                f"File {filename} read successfully, length={len(data)}")
            return data

    except Exception as e:
        logger.error(f"Error reading file {filename}: {e}", exc_info=True)
        raise

def write_file(filename: str, text: str) -> None:
    """
    Writes the given text to a file.

    :param filename: Path to the file to write to.
    :param text: The text
    :return: None
    """

    logger.info(f"Writing results to file: {filename}")

    try:
        with open(filename, 'w', encoding='utf-8') as file:
            file.write(text)
            logger.debug(f"Successfully wrote {len(text)} characters")

    except Exception as e:
        logger.error(f"Error writing file {filename}: {e}", exc_info=True)
        print(f"Error writing file: {e}")


def frequency_bit_test(sequence: str) -> float:
    """
    Performs a frequency bit test

    :param sequence: The bit sequence
    :return: P-value
    """

    logger.debug("Starting frequency bit test")

    n = len(sequence)

    if n == 0:
        logger.warning("Frequency bit test failed: empty sequence")
        raise ValueError("Sequence is empty")

    s = sum([1 if bit == "1" else -1 for bit in sequence])

    p_value = math.erfc((abs(s) / math.sqrt(n)) / math.sqrt(2))

    logger.debug(f"Frequency bit test completed: p-value={p_value}")

    return p_value


def runs_test(sequence: str) -> float:
    """
    Performs a test for identical consecutive bits.

    :param sequence: The sequence
    :return: P-value
    """

    logger.debug("Starting runs test")

    n = len(sequence)
    p = sequence.count('1') / n

    if abs(p - 0.5) >= 2 / math.sqrt(n):
        logger.warning("Runs test condition failed, returning p-value=0")
        return 0.0

    v_n = 0

    for i in range(n - 1):

        if sequence[i] != sequence[i + 1]:
            v_n += 1

    numerator = abs(v_n - 2 * n * p * (1 - p))
    denominator = 2 * math.sqrt(2 * n) * p * (1 - p)

    p_value = math.erfc(numerator / denominator)
    logger.debug(f"Runs test completed: p-value={p_value}")

    return p_value


def block_run_test(sequence: str) -> float:
    """
    Performs a test for the longest sequence of ones in a block.

    :param sequence: The sequence
    :return: P-value
    """

    logger.debug("Starting block run test")

    n = len(sequence)

    if n < 128:
        logger.error("Block run test failed: sequence length < 128")
        raise ValueError("Minimum 128 bits required")

    N = n // 8

    v = [0, 0, 0, 0]

    for i in range(N):

        block = sequence[i * 8 : (i + 1) * 8]

        max_run = 0
        current_run = 0

        for bit in block:

            if bit == '1':
                current_run += 1

                if current_run > max_run:
                    max_run = current_run

            else:
                current_run = 0

        match max_run:

            case 0 | 1:
                v[0] += 1

            case 2:
                v[1] += 1

            case 3:
                v[2] += 1

            case _:
                v[3] += 1

    x_2 = 0.0

    for i in range(len(v)):
        x_2 += (v[i] - 16 * PI[i]) ** 2 / (16 * PI[i])

    p_value = gammainc(3 / 2, x_2 / 2)

    logger.debug(f"Block run test completed: p-value={p_value}")

    return p_value


def main():
    logger.info("Application started")

    cpp_sequence = read_file(cpp_sequence_txt)
    java_sequence = read_file(java_sequence_txt)

    logger.info("Running tests for CPP sequence")

    p_val_freq_bits_cpp = frequency_bit_test(cpp_sequence)
    p_val_ident_bits_cpp = runs_test(cpp_sequence)
    p_val_longest_bits_block_cpp = block_run_test(cpp_sequence)

    result_cpp_test = (f"CPP sequence: {cpp_sequence}\n\n"
                      f"Frequency bit test: {p_val_freq_bits_cpp}\n"
                      f"Test for identical consecutive bits: {p_val_ident_bits_cpp}\n"
                      f"Test for the longest sequence of ones in a block: "
                      f"{p_val_longest_bits_block_cpp}")

    write_file(test_results_cpp, result_cpp_test)

    logger.info("Running tests for Java sequence")

    p_val_freq_bits_java = frequency_bit_test(java_sequence)
    p_val_ident_bits_java = runs_test(java_sequence)
    p_val_longest_bits_block_java = block_run_test(java_sequence)

    result_java_test = (f"Java sequence: {java_sequence}\n\n"
                      f"Frequency bit test: {p_val_freq_bits_java}\n"
                      f"Test for identical consecutive bits: {p_val_ident_bits_java}\n"
                      f"Test for the longest sequence of ones in a block: "
                      f"{p_val_longest_bits_block_java}")

    write_file(test_results_java, result_java_test)

    logger.info("Application finished successfully")


if __name__ == "__main__":
    main()