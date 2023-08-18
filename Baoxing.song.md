
# Run the MCScanX pipeline
The genome sequence and genome annotation data were prepared by Zhiyu Liu
## Extract the protein sequences from the genome sequences and genome annotations for each species
```
gffread -g /media/songlab/songlab_2023_14t/team2/data/Zm-B73-REFERENCE-NAM-5.0.fa -y zea.p.fa /media/songlab/songlab_2023_14t/team2/data/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3
gffread -g /media/songlab/songlab_2023_14t/team2/data/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa -y os.p.fa /media/songlab/songlab_2023_14t/team2/data/Oryza_sativa.IRGSP-1.0.57.gff3
gffread -g /media/songlab/songlab_2023_14t/team2/data/Setaria_viridis.Setaria_viridis_v2.0.dna.toplevel.fa -y sv.p.fa /media/songlab/songlab_2023_14t/team2/data/Setaria_viridis.Setaria_viridis_v2.0.57.gff3
gffread -g /media/songlab/songlab_2023_14t/team2/data/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa -y sb.p.fa /media/songlab/songlab_2023_14t/team2/data/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.57.gff3
```

## Identify and extract the longest pep sequences for each gene
```
python3 achorwave_quota/longestPep.py /media/songlab/songlab_2023_14t/team2/data/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 /media/songlab/songlab_2023_14t/team2/data/Zm-B73-REFERENCE-NAM-5.0.fa zea.p.fa maize.protein.fa
python3 achorwave_quota/longestPep.py /media/songlab/songlab_2023_14t/team2/data/Oryza_sativa.IRGSP-1.0.57.gff3 /media/songlab/songlab_2023_14t/team2/data/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa os.p.fa rice.protein.fa
python3 achorwave_quota/longestPep.py /media/songlab/songlab_2023_14t/team2/data/Setaria_viridis.Setaria_viridis_v2.0.57.gff3 /media/songlab/songlab_2023_14t/team2/data/Setaria_viridis.Setaria_viridis_v2.0.dna.toplevel.fa sv.p.fa setaria.protein.fa
python3 achorwave_quota/longestPep.py /media/songlab/songlab_2023_14t/team2/data/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.57.gff3 /media/songlab/songlab_2023_14t/team2/data/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa sb.p.fa sorghum.protein.fa
```
This python scripts is released at https://github.com/baoxingsong/quota_Anchor


## Perform pep sequences alignment
```
mkdir xyz
/media/songlab/songlab_2023_14t/team2/baoxingsong/blast-2.2.26/bin/formatdb -p T -i maize.protein.fa -n maize.protein
/media/songlab/songlab_2023_14t/team2/baoxingsong/blast-2.2.26/bin/formatdb -p T -i rice.protein.fa -n rice.protein
/media/songlab/songlab_2023_14t/team2/baoxingsong/blast-2.2.26/bin/formatdb -p T -i setaria.protein.fa -n setaria.protein

/media/songlab/songlab_2023_14t/team2/baoxingsong/blast-2.2.26/bin/blastall -i rice.protein.fa -d maize.protein -p blastp -e 1e-10 -b 5 -v 5 -m 8 -o rice.maize.blastp &
/media/songlab/songlab_2023_14t/team2/baoxingsong/blast-2.2.26/bin/blastall -i setaria.protein.fa -d maize.protein -p blastp -e 1e-10 -b 5 -v 5 -m 8 -o setaria.maize.blastp &
/media/songlab/songlab_2023_14t/team2/baoxingsong/blast-2.2.26/bin/blastall -i sorghum.protein.fa -d maize.protein -p blastp -e 1e-10 -b 5 -v 5 -m 8 -o sorghum.maize.blastp &
/media/songlab/songlab_2023_14t/team2/baoxingsong/blast-2.2.26/bin/blastall -i setaria.protein.fa -d rice.protein -p blastp -e 1e-10 -b 5 -v 5 -m 8 -o setaria.rice.blastp &
/media/songlab/songlab_2023_14t/team2/baoxingsong/blast-2.2.26/bin/blastall -i sorghum.protein.fa -d rice.protein -p blastp -e 1e-10 -b 5 -v 5 -m 8 -o sorghum.rice.blastp &
/media/songlab/songlab_2023_14t/team2/baoxingsong/blast-2.2.26/bin/blastall -i sorghum.protein.fa -d setaria.protein -p blastp -e 1e-10 -b 5 -v 5 -m 8 -o sorghum.setaria.blastp &

cat rice.maize.blastp setaria.maize.blastp sorghum.maize.blastp setaria.rice.blastp sorghum.rice.blastp sorghum.setaria.blastp > ./xyz/xyz.blast
```

## Prepare the gene position information for MCScanX
```
cat /media/songlab/songlab_2023_14t/team2/data/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 | grep -P "\tgene\t" | sed '~s/^chr//g' | sed '~s/;bio.*$//g' | sed '~s/ID=//g' | sed '~s/;.*$//g' | awk '{print "Zm"$1"\t"$9"\t"$4"\t"$5}' > ./xyz/xyz.bed
cat /media/songlab/songlab_2023_14t/team2/data/Oryza_sativa.IRGSP-1.0.57.gff3 | grep -P "\tgene\t" | sed '~s/^chr//g' | sed '~s/;bio.*$//g' | sed '~s/ID=//g' | sed '~s/;.*$//g' | awk '{print "Os"$1"\t"$9"\t"$4"\t"$5}' >> ./xyz/xyz.bed
cat /media/songlab/songlab_2023_14t/team2/data/Setaria_viridis.Setaria_viridis_v2.0.57.gff3 | grep -P "\tgene\t" | sed '~s/^chr//g' | sed '~s/;bio.*$//g' | sed '~s/ID=//g' | sed '~s/;.*$//g' | awk '{print "Sv"$1"\t"$9"\t"$4"\t"$5}' >> ./xyz/xyz.bed
cat /media/songlab/songlab_2023_14t/team2/data/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.57.gff3 | grep -P "\tgene\t" | sed '~s/^chr//g' | sed '~s/;bio.*$//g' | sed '~s/ID=//g' | sed '~s/;.*$//g' | awk '{print "Sb"$1"\t"$9"\t"$4"\t"$5}' >> ./xyz/xyz.bed
```

## Run the MCScanX pipeline
```
~/mc/MCScanX/MCScanX ./xyz/xyz
```


# Implement a strand and whole-genome duplication aware syntenic gene identification pipeline
https://github.com/baoxingsong/quota_Anchor

# Attendee
Baoxing Song, baoxing.song@pku-iaas.edu.cn
