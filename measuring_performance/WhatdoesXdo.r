# Checking that cutadapt option X does what we think it does

# Create a probe to find

probe = "TGGAACAGTA"

## Create scenarios

## Shorten probes

shorten_probe <- function (string, number){ str_trunc(string, width = str_length(string) - number, side = "left", ellipsis = "")}

## Move probes right

move_probes <- function(string, number){paste0(paste(rep("A", number),collapse =""), string)}

move_probes(probe, 4)

## AMPLICON

amplicon = "AACGAACGCTGGCGGCAGGCTTAACACATGCAAGTCGAGCGGGCCTCGCAATAGGTGAGCGGCAGACGGGTGAGTAACGCGTGGGAATCATACCTTTTGGTTCGGAACAACTCCGGGAAACTTGAGCTAATACCGGATACGCCCTTACGGGGGAAAGATTTATCGGGGAAAGATGGGCCCGCGTCTGATTAGCTAGTTGGTGGGGTAATGGCTCACCAAGGCGACGATCCGTAGCTGGTCTGAGAGGATGATCAGCCACATTGGGACTGAGACACGGCCCAAACTCCTACGGGAGGCAGCAGTGGGGAATATTGGACAATGGGCGCAAGCCTGATCCAGCCATGCCGCGTGAGTGATGAAGGCCCTAGGGTTGTAAAGCTCTTTTGTGCGGTGAAGATAATGACGGTAACCGCAGAATAAGCCCCGGCTAACTTCGTGCCAGCAGCCGCGGTAATACGAAGGGGGCTAGCGTTGTTCGGAATCACTGGGCGTAAAGGGCGCGTAGGCGGGTCTTTAAGTCAGGGGTGAAATCCTGGAGCTCAACTCCAGAACTGCCTTTGATACTGGGTATCTTGAGTATCGGAAGAGGTGAGTGGAACTGCGAGTGTAGAGGTGAAATTCGTAGATATTCGCAGGAACACCAGTGGCGAAGGCGGCTCACTGGTCCGATTACTGACGCTGAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGAATGTTAGCCGTTGGGGGGTTTACTCTCTAGTGGCGCAGCTAACGCATTAAGCATTCCGCCTGGGGAGTACGGTCGCAAGATTAAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGACGCAACGCGCAGAACCTTACCAGCTCTTGACATGTCGGGGACCGGTCGCAGAGATGTGACCTTCTCCTTCGGTTAGCCTGGAACACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTCGTCCTTAGTTGCCACCATTTAGTTGGGCACTCTAAGGGGACTGCCGGTGATAAGCCGAGAGGAAGGTGGGGATGACGTCAAGTCCTCATGGCCCTTACGGGCTGGGCTACACACGTGCTACAATGGCGGTGACAGTGGGCAGCGAGACAGCGATGTCGAGCTAATCTCCA"

### COMBINE them

tibble(Scenario = c("Shorten", "Moved")) 

Shorten.probes <- tibble (number = 1:4) %>% mutate(Initial.string = map_chr(number,~ shorten_probe(probe, .x))) %>% 
  mutate(cutadapt.input = paste0(Initial.string, amplicon)) %>% mutate(Scenario = "Shorten")

Moved.probes <- tibble(number = 1:20) %>% mutate(Initial.string = map_chr(number, ~ move_probes(probe, .x))) %>% 
  mutate(cutadapt.input = paste0(Initial.string, amplicon)) %>% mutate(Scenario = "Moved")

bind_rows(Shorten.probes, Moved.probes) %>% 
  unite(Scenario, number, col = "Name", sep = "_chars=") -> combined

seqinr::write.fasta(sequences = as.list(combined$cutadapt.input),
                    names = as.list(combined$Name),
                    file.out = "~/cutadapt.test.fasta")
## Run this

system('cutadapt -g "XTGGAACAGTA" -e 0.1 -o ~/cutadapt.output.fasta ~/cutadapt.test.fasta --discard-untrimmed')

# Import result

winners <- readFASTA("~/cutadapt.output.fasta")
winners

# only one of the moved chars, and I think it is because of the -e 0.1

system('cutadapt -g "XTGGAACAGTA" -e 0 -o ~/cutadapt.output.noe.fasta ~/cutadapt.test.fasta --discard-untrimmed')

readFASTA("~/cutadapt.output.noe.fasta")

# YEP
