import re
import cgat.MEME


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
        
        if "E" in self.properties:
            self.evalue = self.properties["E"]
            
        if "P" in self.properties:
            self.evalue = self.properties["P"]

        if "nsites" in self.properties:
            self.nsites = self.properties["nsites"]

        if "w" in self.properties:
            self.motif_len = self.properties["w"]
             
        self.letter_probability_line = line
        

class HomerMotifFile(cgat.MEME.MemeMotifFile):
    matrix_line_re = re.compile(
        "letter-probability matrix: .+w=\s?([0-9]+).+[EP]=\s?(.+)")
