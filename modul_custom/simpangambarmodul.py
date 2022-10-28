from rdkit.Chem import Draw
def simpangambar(simpan, catatan):
    varsimpan = Draw.MolToImage(simpan)
    
    nourut = open("D:/dokumen/Academic/Pemograman/python/Library RDkit/RDkit output/nourut.txt", "r")
    lines = nourut.readlines()
    nolines = lines[0]
    splitno = nolines.split(',')
    banyakno = len(splitno)

    tulis = banyakno + 1
    strtulis = str(tulis)+ ","
    direktorisimpan = f"D:/dokumen/Academic/Pemograman/python/Library RDkit/RDkit output/nourut.txt"
    tulisnourut = open(direktorisimpan, "a")
    tulis = tulisnourut.write(strtulis)
    pesan = print(f"gambar disimpan didekat: {direktorisimpan}")

    varsimpan.save(f"D:/dokumen/Academic/Pemograman/python/Library RDkit/RDkit output/gambar{strtulis} - {catatan}.png")
    return pesan