from Bio import motifs 



with open("/home/nadzhou/SEQs/navid/meme12/meme.xml", "r") as handle: 
    m = motifs.parse(handle, "meme")


print(m)


for motif in m: 
    print(motif)