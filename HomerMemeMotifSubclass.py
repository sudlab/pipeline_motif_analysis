import re
import cgat.MEME
from cgat.MEME import MotifList


class HomerMotif(cgat.MEME.MemeMotif):

    def parse_probability_line(self, line):

        line = line.strip()
        properties = re.findall("(\S+)=\s?(\S+)", line)

        for key, value in properties:
            print("key pair is %s, %s" % (key,value))
            try:
                value = int(value)
            except ValueError:
                try:
                    value = float(value)
                except ValueError:
                    pass

            self.properties[key] = value

        if "P" in self.properties:
            self.evalue = self.properties["P"]

        if "E" in self.properties:
            self.evalue = self.properties["E"]

        if "nsites" in self.properties:
            self.nsites = self.properties["nsites"]

        if "w" in self.properties:
            self.motif_len = self.properties["w"]

        self.letter_probability_line = line


class HomerMotifFile(cgat.MEME.MemeMotifFile):
    matrix_line_re = re.compile(
        "letter-probability matrix: .+w=\s?([0-9]+).+[EP]=\s?(.+)")

    def __init__(self, buffer_or_Motiffile):

        MotifList.__init__(self)

        try:
            self.header = buffer_or_Motiffile.header
        except AttributeError:

            self.header = []
            current_motif = None
            for line in buffer_or_Motiffile:
                line.strip()
                if self.motif_line_re.match(line):
                    current_motif = HomerMotif(line)
                    break
                else:
                    self.header.append(line)

            self.header = '\n'.join(self.header)

            if current_motif is None:
                return

            for line in buffer_or_Motiffile:
                line.strip()

                if self.motif_line_re.match(line):
                    self.motifs.append(current_motif)
                    current_motif = HomerMotif(line)

                if self.matrix_line_re.match(line):
                    current_motif.parse_probability_line(line)
                    for i in range(current_motif.motif_len):
                        line = next(buffer_or_Motiffile)
                        line.strip()
                        current_motif.letter_probability_matrix.append(
                            line.split())
            else:
                self.motifs.append(current_motif)
