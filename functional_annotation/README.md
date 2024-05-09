# 0. Building SnpEff database

The Iberian lynx reference genome is not found in the SnpEff database. We will build the database using the reference genome and the annotation file. 
For that, we first need to modify the snpEff.config file to include the new database. As I don't have the permissions to write the file, I need to create a copy in my $STORE.


```bash
module load snpeff/5.0
cp /opt/cesga/2020/software/Core/snpeff/5.0/snpEff.config /mnt/netapp1/Store_CSIC/home/csic/eye/lmf/snpEff/
```

I need to add these lines at the end of the config file:
```bash
# Lynx_pardinus, version mLynPar1.2
LYPA1_2A.genome : Iberian lynx        # from now on, LYPA1_2A is the code for the Lynx pardinus reference genome (in snpEff)
```
