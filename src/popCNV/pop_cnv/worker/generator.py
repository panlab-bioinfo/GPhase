class BIN:
    def __init__(self):
        pass

    @staticmethod
    def generate(bin_size, genome_db):
        bins = {}
        for chrn in genome_db:
            bins[chrn] = []
            i = 1
            while i <= len(genome_db[chrn]):
                e = i + bin_size - 1
                if e >= len(genome_db[chrn]):
                    e = len(genome_db[chrn])
                bins[chrn].append([i, e])
                i += bin_size
        return bins
