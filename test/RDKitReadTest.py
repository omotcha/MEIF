"""
platform: win
env: any with RDKit
name: RDKitReadTest.py
RDKit read molecule Tester
"""
import os

from rdkit import Chem
from configs.config import dataset_2016_dir, casf_dir
from util.RDKitHelper import Mol2MolSupplier_test, get_decoy_names

block_good = """
@<TRIPOS>MOLECULE
1c5z_415
   18    18     1     1     0
SMALL
MMFF94_CHARGES


@<TRIPOS>ATOM
      1 C1          9.1430    0.8640   23.8760 C.ar      1 BAM        -0.0220 
      2 C2          8.6700   -0.1080   24.7830 C.ar      1 BAM        -0.1500 
      3 C3          8.3120    0.2320   26.0860 C.ar      1 BAM        -0.1500 
      4 C4          8.4170    1.5490   26.5100 C.ar      1 BAM        -0.1500 
      5 C5          8.8790    2.5250   25.6470 C.ar      1 BAM        -0.1500 
      6 C6          9.2330    2.1830   24.3550 C.ar      1 BAM        -0.1500 
      7 C           9.5040    0.5340   22.5550 C.cat     1 BAM         0.7308 
      8 N1          8.7100   -0.2060   21.8090 N.pl3     1 BAM        -1.1924 
      9 N2         10.6390    0.9750   22.0840 N.pl3     1 BAM        -1.1924 
     10 H1          8.5830   -1.1390   24.4600 H         1 BAM         0.1500 
     11 H2          7.9520   -0.5310   26.7660 H         1 BAM         0.1500 
     12 H3          8.1360    1.8150   27.5220 H         1 BAM         0.1500 
     13 H4          8.9630    3.5530   25.9810 H         1 BAM         0.1500 
     14 H5          9.5930    2.9580   23.6890 H         1 BAM         0.1500 
     15 H6          8.9830   -0.4430   20.8440 H         1 BAM         0.6690 
     16 H7          7.8170   -0.5490   22.1900 H         1 BAM         0.6690 
     17 H8         10.9220    0.7450   21.1210 H         1 BAM         0.6690 
     18 H9         11.2530    1.5540   22.6770 H         1 BAM         0.6690 
@<TRIPOS>BOND
     1    1    2 ar   
     2    1    6 ar   
     3    1    7 1    
     4    2    3 ar   
     5    2   10 1    
     6    3    4 ar   
     7    3   11 1    
     8    4    5 ar   
     9    4   12 1    
    10    5    6 ar   
    11    5   13 1    
    12    6   14 1    
    13    7    8 ar   
    14    7    9 ar   
    15    8   15 1    
    16    8   16 1    
    17    9   17 1    
    18    9   18 1    
@<TRIPOS>SUBSTRUCTURE
     1 BAM         1 GROUP             4 ****  ****    0 ROOT 
@<TRIPOS>NORMAL
@<TRIPOS>ALT_TYPE
MMFF94_ALT_TYPE_SET 
MMFF94 1 CB 2 CB 3 CB 4 CB 5 CB 6 CB 7 CNN+ 10 HC 11 HC 12 HC 13 HC \
14 HC 15 HN 16 HN 17 HN 18 HN 8 NCN+ 9 NCN+ 
"""

block_bad = """
@<TRIPOS>MOLECULE
1c5z_247
   18    18     1     1     2
SMALL
MMFF94_CHARGES


@<TRIPOS>ATOM
      1 C1          8.7877    0.6795   25.0176 C.ar      1 BAM        -0.0220 
      2 C2          7.9780    1.7432   25.4686 C.ar      1 BAM        -0.1500 
      3 C3          8.0583    3.0089   24.8917 C.ar      1 BAM        -0.1500 
      4 C4          8.9468    3.2391   23.8510 C.ar      1 BAM        -0.1500 
      5 C5          9.7555    2.2190   23.3868 C.ar      1 BAM        -0.1500 
      6 C6          9.6729    0.9652   23.9622 C.ar      1 BAM        -0.1500 
      7 C           8.7125   -0.6049   25.5887 C.cat      1 BAM         0.7308 
      8 N1          8.6322   -1.6741   24.8251 N.pl3     1 BAM        -1.1924 
      9 N2          8.7236   -0.7156   26.8904 N.pl3     1 BAM        -1.1924 
     10 H1          7.2798    1.5725   26.2800 H         1 BAM         0.1500 
     11 H2          7.4279    3.8123   25.2551 H         1 BAM         0.1500 
     12 H3          9.0069    4.2232   23.4002 H         1 BAM         0.1500 
     13 H4         10.4509    2.4013   22.5755 H         1 BAM         0.1500 
     14 H5         10.3124    0.1743   23.5873 H         1 BAM         0.1500 
     15 H6          8.5776   -2.6096   25.2528 H         1 BAM         0.6690 
     16 H7          8.6240   -1.5753   23.7999 H         1 BAM         0.6690 
     17 H8          8.7868    0.1276   27.4785 H         1 BAM         0.6690 
     18 H9          8.6692   -1.6462   27.3287 H         1 BAM         0.6690 
@<TRIPOS>BOND
     1    1    2 ar   
     2    1    6 ar   
     3    1    7 1    
     4    2    3 ar   
     5    3    4 ar   
     6    4    5 ar   
     7    5    6 ar   
     8    7    8 ar   
     9    7    9 ar   
    10    2   10 1    
    11    3   11 1    
    12    4   12 1    
    13    5   13 1    
    14    6   14 1    
    15    8   15 1    
    16    8   16 1    
    17    9   17 1    
    18    9   18 1    
@<TRIPOS>SUBSTRUCTURE
     1 BAM         1 PERM              0 ****  ****    0 ROOT 
@<TRIPOS>SET
DONOR_HYDROGENS STATIC     ATOMS    <user>   **** ""
4 15 16 17 18
ATOM$RED        STATIC     ATOMS    COLORGROUP SYSTEM 
4 15 16 17 18
@<TRIPOS>NORMAL
@<TRIPOS>ALT_TYPE
MMFF94_ALT_TYPE_SET 
MMFF94 1 CB 2 CB 3 CB 4 CB 5 CB 6 CB 7 CNN+ 10 HC 11 HC 12 HC 13 HC \
14 HC 15 HN 16 HN 17 HN 18 HN 8 NCN+ 9 NCN+ 
"""

block_test = """
@<TRIPOS>MOLECULE
1c5z_247
   18    18     1     1     2
SMALL
MMFF94_CHARGES


@<TRIPOS>ATOM
      1 C1          8.7877    0.6795   25.0176 C.ar      1 BAM        -0.0220 
      2 C2          7.9780    1.7432   25.4686 C.ar      1 BAM        -0.1500 
      3 C3          8.0583    3.0089   24.8917 C.ar      1 BAM        -0.1500 
      4 C4          8.9468    3.2391   23.8510 C.ar      1 BAM        -0.1500 
      5 C5          9.7555    2.2190   23.3868 C.ar      1 BAM        -0.1500 
      6 C6          9.6729    0.9652   23.9622 C.ar      1 BAM        -0.1500 
      7 C           8.7125   -0.6049   25.5887 C.ar      1 BAM         0.7308 
      8 N1          8.6322   -1.6741   24.8251 N.pl3     1 BAM        -1.1924 
      9 N2          8.7236   -0.7156   26.8904 N.pl3     1 BAM        -1.1924 
     10 H1          7.2798    1.5725   26.2800 H         1 BAM         0.1500 
     11 H2          7.4279    3.8123   25.2551 H         1 BAM         0.1500 
     12 H3          9.0069    4.2232   23.4002 H         1 BAM         0.1500 
     13 H4         10.4509    2.4013   22.5755 H         1 BAM         0.1500 
     14 H5         10.3124    0.1743   23.5873 H         1 BAM         0.1500 
     15 H6          8.5776   -2.6096   25.2528 H         1 BAM         0.6690 
     16 H7          8.6240   -1.5753   23.7999 H         1 BAM         0.6690 
     17 H8          8.7868    0.1276   27.4785 H         1 BAM         0.6690 
     18 H9          8.6692   -1.6462   27.3287 H         1 BAM         0.6690 
@<TRIPOS>BOND
     1    1    2 ar   
     2    1    6 ar   
     3    1    7 1    
     4    2    3 ar   
     5    3    4 ar   
     6    4    5 ar   
     7    5    6 ar   
     8    7    8 ar   
     9    7    9 ar   
    10    2   10 1    
    11    3   11 1    
    12    4   12 1    
    13    5   13 1    
    14    6   14 1    
    15    8   15 1    
    16    8   16 1    
    17    9   17 1    
    18    9   18 1    
@<TRIPOS>SUBSTRUCTURE
     1 BAM         1 PERM              0 ****  ****    0 ROOT 
@<TRIPOS>SET
DONOR_HYDROGENS STATIC     ATOMS    <user>   **** ""
4 15 16 17 18
ATOM$RED        STATIC     ATOMS    COLORGROUP SYSTEM 
4 15 16 17 18
@<TRIPOS>NORMAL
@<TRIPOS>ALT_TYPE
MMFF94_ALT_TYPE_SET 
MMFF94 1 CB 2 CB 3 CB 4 CB 5 CB 6 CB 7 CNN+ 10 HC 11 HC 12 HC 13 HC \
14 HC 15 HN 16 HN 17 HN 18 HN 8 NCN+ 9 NCN+ 
"""


def testReadFromBlock():
    m = Chem.MolFromMol2Block(block_bad, sanitize=True)
    print(m.GetNumAtoms())


def testReadFromMol2():
    lids = os.listdir(dataset_2016_dir["core"])
    bad_lid = []
    for lid in lids:
        m = Chem.MolFromMol2File(os.path.join(dataset_2016_dir["core"], lid, "{}_ligand.mol2".format(lid)))
        if m is None:
            bad_lid.append(lid)
    print(bad_lid)


def testBadNum(f):
    mols = Mol2MolSupplier_test(f, sanitize=True)
    names = get_decoy_names(f)
    print("{}: {} can be identified from {} decoys".format(f[-16:-12], len(mols), len(names)))


if __name__ == '__main__':
    # testBadNum(os.path.join(casf_dir["decoys"], "1c5z_decoys.mol2"))
    testReadFromBlock()

