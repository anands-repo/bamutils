import argparse
import os
import logging
from split_bams import Categorical

def close_file(f):
    if not f.closed:
        f.close()


def read_single_fastq(fq):
    text = ""
    for i in range(4):
        text += next(fq)
    return text


class FastqReader:
    def __init__(self, fq1, fq2=None):
        self.fq1 = fq1
        self.fq2 = fq2

    def _close(self):
        if hasattr(self, "_fhandle0"):
            close_file(self._fhandle0)

        if hasattr(self, '_fhandle1'):
            close_file(self._fhandle1)

    def __iter__(self):
        self._close()

        self._fhandle0 = open(self.fq1, 'r')

        if self.fq2:
            self._fhandle1 = open(self.fq2, 'r')

        return self

    def __next__(self):
        lines = [read_single_fastq(self._fhandle0)]
        if hasattr(self, '_fhandle1'):
            lines.append(read_single_fastq(self._fhandle1))
        return lines


class FastqWriter:
    def __init__(self, prefix, paired_end=False):
        self.prefix = prefix

        if os.path.exists(self.prefix):
            raise ValueError("Directory %s exists" % self.prefix)

        os.makedirs(self.prefix)
        self.paired_end = paired_end
        self._fhandle0 = open(os.path.join(self.prefix, "1.fq"), "w")

        if self.paired_end:
            self._fhandle1 = open(os.path.join(self.prefix, "2.fq"), "w")

    def write(self, lines):
        self._fhandle0.write(lines[0])

        if len(lines) == 2:
            assert(self.paired_end)
            self._fhandle1.write(lines[1])

        if self.paired_end:
            assert(len(lines) == 2)

    def close(self):
        if hasattr(self, '_fhandle0'):
            close_file(self._fhandle0)

        if hasattr(self, '_fhandle1'):
            close_file(self._fhandle1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Split a source FASTQ file into multiple parts")

    parser.add_argument(
        "--source",
        help="Source fastq files (comma-separated for paired-end)",
        required=True,
    )

    parser.add_argument(
        "--prefix",
        help="Prefix of output files",
        required=True,
    )

    parser.add_argument(
        "--proportions",
        help="Comma-separated proportions of lines to be sent to each file",
        required=True,
    )

    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format='%(asctime)-15s %(message)s')
    proportions = list(map(float, args.proportions.split(",")))
    distribution = Categorical(proportions)

    source = args.source.split(",")
    reader = FastqReader(*source)
    output_prefixes = [
        args.prefix + "_%d" % i for i in range(len(proportions))
    ]
    ohandles = {
        op: FastqWriter(op, len(source) == 2) for op in output_prefixes
    }
    num_entries = {
        op: 0 for op in output_prefixes
    }

    for i, lines in enumerate(iter(reader)):
        target = output_prefixes[distribution.sample()]
        thandle = ohandles[target]
        num_entries[target] += 1
        thandle.write(lines)

        if (i + 1) % 1000000 == 0:
            logging.info("Completed %d reads" % (i + 1))
            logging.info("Entries per file = %s" % str(num_entries))

    for o in ohandles.values():
        o.close()
