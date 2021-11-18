#!/bin/bash

#SBATCH -A snic2021-5-28
#SBATCH -J RepeatModels
#SBATCH -o RepeatModels_%j.out
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 1-00:00:00

### Description
# A pipeine developped specifically to build a substitution model based on repeats calculated on the most ancestral sequence of the 241 Mammals Alignment.
# Designed to run on the Uppmax	SLURM system. Can be simply edited.
# Script customized to run everything in one run.
# if not avaliable as module on your system, the following packages/ apps are required to be installed :
# You will need to change the command lines in accordance.
# - Repeatmasker : https://www.repeatmasker.org/RepeatMasker/
# - hal : https://github.com/ComparativeGenomicsToolkit/hal
# - BEDOPS v2.4.38 : https://bedops.readthedocs.io/en/latest
# - BEDTools v2.29.2 : https://bedtools.readthedocs.io/en/latest
# - Phast v1.5 : https://github.com/CshlSiepelLab/phast
# - mafTools from Dent Earl : https://github.com/dentearl/mafTools

echo -n "Time started: "
date

# Load modules available on Uppmax
module load bioinfo-tools
module load hal
module load RepeatMasker
module load BEDOPS
module load BEDTools
module load phast/v1.5


### Extract Anc110 sequence
hal_file="/proj/uppstore2017228/KLT.04.200M/200m_MD/data/new_250_MAMMALS_v2_20201120/HAL/241-mammalian-2020v2.hal"
species_name="fullTreeAnc110"
out_dir="/proj/uppstore2017228/KLT.04.200M/200m_MD/data/new_250_MAMMALS_v2_20201120/HAL"
out_file=$out_dir"/"$species_name".fa"
hal2fasta $hal_file $species_name > $out_file

### RepeatMask the ancestral sequence of reference
RepeatMasker $out_file # Output will be in folder with the same name than ref.

### Filter out problematic repeats # Credits to A.Smit
rmskOUT=$out_dir"/"$species_name"/$species_name".fa.out"
rmskOUTF=$out_dir"/"$species_name"/$species_name".filtered.fa.out"
rmskOUTBEDF=$out_dir"/"$species_name"/$species_name".filtered.merged.fa.bed"
./exclude.pl $rmskOUT > $rmskOUTF
awk -v OFS='\t' '{print $5,$6,$7}' $rmskOUTF | sort -k1,1 -k2,2n | bedtools merge -i - > $rmskOUTBEDF

### Check synteny with some branches of the Primates tree:
infile=$out_dir"/LO_"$species_name".bed"
final_file=$out_dir"/$species_name"_Allfam.bed"
list_ancestor="fullTreeAnc112 fullTreeAnc107 fullTreeAnc103 fullTreeAnc88 fullTreeAnc78 fullTreeAnc70"
for ancestral_branch in $list_ancestor ; do
	outfile=$out_dir"/.filtered.merged.vs"$ancestral_branch".fa.bed"
	halLiftover $halfile $species2 $infile $ancestral_branch $outfile --noDupes --outPSL
	awk -v OFS='\t' '{print $10,$12,$13,$14,$16,$17}' $outfile > $outfile"_temp" && mv $outfile"_temp" $outfile
	bedtools intersect -a $final_file -b $outfile > $final_file"_temp" && mv $final_file"_temp" $final_file
done

### Calculate models

# Calculate model on random positions of ancestral repeat coordinates
out_AR100kb=$out_dir"/random_100kb_AncRepeats.bed"
bedops --chop 1 $final_file | shuf -n 100000 | sort -k1,1 -k2,2n | bedtools merge -i - > $out_AR100kb

#hal2maf
maf_100kb=$out_dir"/random_100kb_AncRepeats.maf"
maf_100kbmod=$out_dir"/random_100kb_AncRepeats.mod"
genomes="Acinonyx_jubatus,Acomys_cahirinus,Ailuropoda_melanoleuca,Ailurus_fulgens,Allactaga_bullata,Alouatta_palliata,Ammotragus_lervia,Anoura_caudifer,Antilocapra_americana,Aotus_nancymaae,Aplodontia_rufa,Artibeus_jamaicensis,Ateles_geoffroyi,Balaenoptera_acutorostrata,Balaenoptera_bonaerensis,Beatragus_hunteri,Bison_bison,Bos_indicus,Bos_mutus,Bos_taurus,Bubalus_bubalis,Callicebus_donacophilus,Callithrix_jacchus,Camelus_bactrianus,Camelus_dromedarius,Camelus_ferus,Canis_lupus,Canis_lupus_familiaris,Capra_aegagrus,Capra_hircus,Capromys_pilorides,Carollia_perspicillata,Castor_canadensis,Catagonus_wagneri,Cavia_aperea,Cavia_porcellus,Cavia_tschudii,Cebus_albifrons,Cebus_capucinus,Ceratotherium_simum,Ceratotherium_simum_cottoni,Cercocebus_atys,Cercopithecus_neglectus,Chaetophractus_vellerosus,Cheirogaleus_medius,Chinchilla_lanigera,Chlorocebus_sabaeus,Choloepus_didactylus,Choloepus_hoffmanni,Chrysochloris_asiatica,Colobus_angolensis,Condylura_cristata,Craseonycteris_thonglongyai,Cricetomys_gambianus,Cricetulus_griseus,Crocidura_indochinensis,Cryptoprocta_ferox,Ctenodactylus_gundi,Ctenomys_sociabilis,Cuniculus_paca,Dasyprocta_punctata,Dasypus_novemcinctus,Daubentonia_madagascariensis,Delphinapterus_leucas,Desmodus_rotundus,Dicerorhinus_sumatrensis,Diceros_bicornis,Dinomys_branickii,Dipodomys_ordii,Dipodomys_stephensi,Dolichotis_patagonum,Echinops_telfairi,Eidolon_helvum,Elaphurus_davidianus,Elephantulus_edwardii,Ellobius_lutescens,Ellobius_talpinus,Enhydra_lutris,Eptesicus_fuscus,Equus_asinus,Equus_caballus,Equus_przewalskii,Erinaceus_europaeus,Erythrocebus_patas,Eschrichtius_robustus,Eubalaena_japonica,Eulemur_flavifrons,Eulemur_fulvus,Felis_catus,Felis_nigripes,Fukomys_damarensis,Galeopterus_variegatus,Giraffa_tippelskirchi,Glis_glis,Gorilla_gorilla,Graphiurus_murinus,Helogale_parvula,Hemitragus_hylocrius,Heterocephalus_glaber,Heterohyrax_brucei,Hippopotamus_amphibius,Hipposideros_armiger,Hipposideros_galeritus,Homo_sapiens,Hyaena_hyaena,Hydrochoerus_hydrochaeris,Hystrix_cristata,Ictidomys_tridecemlineatus,Indri_indri,Inia_geoffrensis,Jaculus_jaculus,Kogia_breviceps,Lasiurus_borealis,Lemur_catta,Leptonychotes_weddellii,Lepus_americanus,Lipotes_vexillifer,Loxodonta_africana,Lycaon_pictus,Macaca_fascicularis,Macaca_mulatta,Macaca_nemestrina,Macroglossus_sobrinus,Mandrillus_leucophaeus,Manis_javanica,Manis_pentadactyla,Marmota_marmota,Megaderma_lyra,Mellivora_capensis,Meriones_unguiculatus,Mesocricetus_auratus,Mesoplodon_bidens,Microcebus_murinus,Microgale_talazaci,Micronycteris_hirsuta,Microtus_ochrogaster,Miniopterus_natalensis,Miniopterus_schreibersii,Mirounga_angustirostris,Mirza_coquereli,Monodon_monoceros,Mormoops_blainvillei,Moschus_moschiferus,Mungos_mungo,Murina_feae,Muscardinus_avellanarius,Mus_caroli,Mus_musculus,Mus_pahari,Mus_spretus,Mustela_putorius,Myocastor_coypus,Myotis_brandtii,Myotis_davidii,Myotis_lucifugus,Myotis_myotis,Myrmecophaga_tridactyla,Nannospalax_galili,Nasalis_larvatus,Neomonachus_schauinslandi,Neophocaena_asiaeorientalis,Noctilio_leporinus,Nomascus_leucogenys,Nycticebus_coucang,Ochotona_princeps,Octodon_degus,Odobenus_rosmarus,Odocoileus_virginianus,Okapia_johnstoni,Ondatra_zibethicus,Onychomys_torridus,Orcinus_orca,Orycteropus_afer,Oryctolagus_cuniculus,Otolemur_garnettii,Ovis_aries,Ovis_canadensis,Pan_paniscus,Panthera_onca,Panthera_pardus,Panthera_tigris,Pantholops_hodgsonii,Pan_troglodytes,Papio_anubis,Paradoxurus_hermaphroditus,Perognathus_longimembris,Peromyscus_maniculatus,Petromus_typicus,Phocoena_phocoena,Piliocolobus_tephrosceles,Pipistrellus_pipistrellus,Pithecia_pithecia,Platanista_gangetica,Pongo_abelii,Procavia_capensis,Propithecus_coquereli,Psammomys_obesus,Pteronotus_parnellii,Pteronura_brasiliensis,Pteropus_alecto,Pteropus_vampyrus,Puma_concolor,Pygathrix_nemaeus,Rangifer_tarandus,Rattus_norvegicus,Rhinolophus_sinicus,Rhinopithecus_bieti,Rhinopithecus_roxellana,Rousettus_aegyptiacus,Saguinus_imperator,Saiga_tatarica,Saimiri_boliviensis,Scalopus_aquaticus,Semnopithecus_entellus,Sigmodon_hispidus,Solenodon_paradoxus,Sorex_araneus,Spermophilus_dauricus,Spilogale_gracilis,Suricata_suricatta,Sus_scrofa,Tadarida_brasiliensis,Tamandua_tetradactyla,Tapirus_indicus,Tapirus_terrestris,Thryonomys_swinderianus,Tolypeutes_matacus,Tonatia_saurophila,Tragulus_javanicus,Trichechus_manatus,Tupaia_chinensis,Tupaia_tana,Tursiops_truncatus,Uropsilus_gracilis,Ursus_maritimus,Vicugna_pacos,Vulpes_lagopus,Xerus_inauris,Zalophus_californianus,Zapus_hudsonius,Ziphius_cavirostris,fullTreeAnc239"
options="--refGenome fullTreeAnc239 --onlyOrthologs --targetGenomes $genomes"
hal2maf $options --refTargets $out_AR100kb $halfile $maf_100kb
/proj/uppstore2017228/KLT.04.200M/200m_MD/resources/mafTools/bin/mafDuplicateFilter -m $maf_100kb > $maf_100kb"_temp" && mv $maf_100kb"_temp" $maf_100kb

#phyloFit
tree="/proj/uppstore2017228/KLT.04.200M/200m_MD/data/new_250_MAMMALS_v2_20201120/TREE/241-mammalian-2020v2.nh"
phyloFit $maf_100kb -i MAF --subst-mod REV --EM --tree $tree --out-root $maf_100kbmod

### Chr X and Y models
# For models on chrX and Y, conversion of coordinates to a species with identified chrX and Y 
# In that case, we picked human.
species_reference="Homo_sapiens"
final_file_LO=$out_dir"/AncRepeats_Allfam_LO.bed"
halLiftOver $hal_file $species_name $final_file $species_reference $final_file_LO --outPSL --noDupes

# Select random 100kb on chrX and chrY
out_AR100kb_hg38X=$out_dir"/random_100kb_AncRepeats_hg38X.bed"
out_AR100kb_hg38Y=$out_dir"/random_100kb_AncRepeats_hg38Y.bed"
awk -v OFS='\t' '{if ($1=="chrX") print}' $final_file | bedops --chop 1 - | shuf -n 100000 | sort -k1,1 -k2,2n | bedtools merge -i - > $out_AR100kb_hg38X
awk -v OFS='\t' '{if ($1=="chrY") print}' $final_file | bedops --chop 1 - | shuf -n 100000 | sort -k1,1 -k2,2n | bedtools merge -i - > $out_AR100kb_hg38Y

# Extract / calculate model for X: 
maf_100kb_X=$out_dir"/random_100kb_AncRepeats_hg38X.maf"
maf_100kbmod_X=$out_dir"/random_100kb_AncRepeats_hg38X.mod"
options="--refGenome $species_reference --onlyOrthologs --targetGenomes $genomes"
hal2maf $options --refTargets $out_AR100kb_hg38X $hal_file $maf_100kb_X
phyloFit $maf_100kb_X -i MAF --subst-mod REV --EM --tree $tree --out-root $maf_100kbmod_X

# Extract / calculate model for Y:
maf_100kb_Y=$out_dir"/random_100kb_AncRepeats_hg38Y.maf"
maf_100kbmod_Y=$out_dir"/random_100kb_AncRepeats_hg38Y.mod"
options="--refGenome $species_reference --onlyOrthologs --targetGenomes $genomes"
hal2maf $options --refTargets $out_AR100kb_hg38Y $hal_file $maf_100kb_Y
phyloFit $maf_100kb_Y -i MAF --subst-mod REV --EM --tree $tree --out-root $maf_100kbmod_Y

echo -n "Time ended: "
date
