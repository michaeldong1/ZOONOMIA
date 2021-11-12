#!/bin/bash

#SBATCH -A projectName
#SBATCH -J RepeatModel_%j
#SBATCH -o error_out/RepeatModel_%j.out
#SBATCH -o error_out/RepeatModel_%j.err
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 1-00:00:00

### Description
# A pipeine developped specifically to build a substitution model based on repeats calculated on the most ancestral sequence of the 241 Mammals Alignment.
# Following recoomentations, we calculated repeats on the 2nd most ancestral sequence before reconverting coordinates to most ancestral consensus.
# Designed to run on the Uppmax	SLURM system. Can be simply ed
# Script customized to run everything in one run.
# !!! Run the submit script instead. README in run_modelGeneration_mammals.sh

echo -n "Time started: "
date

# Load modules available on Uppmax
module load bioinfo-tools
module load hal
module load RepeatMasker
module load BEDOPS
module load BEDTools
module load phast/v1.5

### PARAMETERS

# Files
hal_file="HAL/241-mammalian-2020v2.hal" # pathway to HAL file
tree="input/TREE/Zoonomia_ChrX_lessGC40_241species_30Consensus.nh" # Newick format tree

# species anc ancestors
species_rmsk="fullTreeAnc238" # Ancestor to calculate repeats on. Can be different from most ancestral one.
species_final="fullTreeAnc239" # most ancestal species. If identical to "species_rmsk", leave blank.
if [[ "$species_final" != "" ]]; then
	species_final=$species_rmsk
fi
list_ancestor="fullTreeAnc5 fullTreeAnc13 fullTreeAnc115 fullTreeAnc237" # List of ancestors or species to check synteny against.
species_reference="Homo_sapiens" # Species of reference for X and Y models.

# Output directory
output_dir="files/output/model_"$species_final # Output directory (will be generated if it does not exist)
if [ ! -d "$output_dir" ]; then
 	mkdir -p "$output_dir" || exit
fi

### output files
# rmsk outputs
speciesrmsk_fasta=$output_dir"/rmsk."$species_rmsk".fa"
rmskOUT=$output_dir"/rmsk."$species_rmsk".fa.out"
rmskOUTF=$output_dir"/rmsk."$species_rmsk".filtered.fa.out"
rmskOUTBEDF=$output_dir"/rmsk."$species_rmsk".filtered.merged.fa.bed"

# name of output if coordinates need reconversion to final
ancestral_repeat_final=$output_dir"/rmsk.LO_"$species_final".bed"
ancestral_repeat_final_synt=$output_dir"/rmsk.LO_"$species_final".Allfam.bed"
ancestral_repeat_final_AR100kb=$output_dir"/rmsk.LO_"$species_final".Allfam.random100kb.bed" # name of random 100k positions 

# hal2maf parameters : extract positions
maf_100kb=$output_dir"/rmsk.LO_"$species_final".Allfam.random100kb.maf"
maf_100kbmod=$output_dir"/rmsk.LO_"$species_final".Allfam.random100kb.mod"
genomes="Acinonyx_jubatus,Acomys_cahirinus,Ailuropoda_melanoleuca,Ailurus_fulgens,Allactaga_bullata,Alouatta_palliata,Ammotragus_lervia,Anoura_caudifer,Antilocapra_americana,Aotus_nancymaae,Aplodontia_rufa,Artibeus_jamaicensis,Ateles_geoffroyi,Balaenoptera_acutorostrata,Balaenoptera_bonaerensis,Beatragus_hunteri,Bison_bison,Bos_indicus,Bos_mutus,Bos_taurus,Bubalus_bubalis,Callicebus_donacophilus,Callithrix_jacchus,Camelus_bactrianus,Camelus_dromedarius,Camelus_ferus,Canis_lupus,Canis_lupus_familiaris,Capra_aegagrus,Capra_hircus,Capromys_pilorides,Carollia_perspicillata,Castor_canadensis,Catagonus_wagneri,Cavia_aperea,Cavia_porcellus,Cavia_tschudii,Cebus_albifrons,Cebus_capucinus,Ceratotherium_simum,Ceratotherium_simum_cottoni,Cercocebus_atys,Cercopithecus_neglectus,Chaetophractus_vellerosus,Cheirogaleus_medius,Chinchilla_lanigera,Chlorocebus_sabaeus,Choloepus_didactylus,Choloepus_hoffmanni,Chrysochloris_asiatica,Colobus_angolensis,Condylura_cristata,Craseonycteris_thonglongyai,Cricetomys_gambianus,Cricetulus_griseus,Crocidura_indochinensis,Cryptoprocta_ferox,Ctenodactylus_gundi,Ctenomys_sociabilis,Cuniculus_paca,Dasyprocta_punctata,Dasypus_novemcinctus,Daubentonia_madagascariensis,Delphinapterus_leucas,Desmodus_rotundus,Dicerorhinus_sumatrensis,Diceros_bicornis,Dinomys_branickii,Dipodomys_ordii,Dipodomys_stephensi,Dolichotis_patagonum,Echinops_telfairi,Eidolon_helvum,Elaphurus_davidianus,Elephantulus_edwardii,Ellobius_lutescens,Ellobius_talpinus,Enhydra_lutris,Eptesicus_fuscus,Equus_asinus,Equus_caballus,Equus_przewalskii,Erinaceus_europaeus,Erythrocebus_patas,Eschrichtius_robustus,Eubalaena_japonica,Eulemur_flavifrons,Eulemur_fulvus,Felis_catus,Felis_nigripes,Fukomys_damarensis,Galeopterus_variegatus,Giraffa_tippelskirchi,Glis_glis,Gorilla_gorilla,Graphiurus_murinus,Helogale_parvula,Hemitragus_hylocrius,Heterocephalus_glaber,Heterohyrax_brucei,Hippopotamus_amphibius,Hipposideros_armiger,Hipposideros_galeritus,Homo_sapiens,Hyaena_hyaena,Hydrochoerus_hydrochaeris,Hystrix_cristata,Ictidomys_tridecemlineatus,Indri_indri,Inia_geoffrensis,Jaculus_jaculus,Kogia_breviceps,Lasiurus_borealis,Lemur_catta,Leptonychotes_weddellii,Lepus_americanus,Lipotes_vexillifer,Loxodonta_africana,Lycaon_pictus,Macaca_fascicularis,Macaca_mulatta,Macaca_nemestrina,Macroglossus_sobrinus,Mandrillus_leucophaeus,Manis_javanica,Manis_pentadactyla,Marmota_marmota,Megaderma_lyra,Mellivora_capensis,Meriones_unguiculatus,Mesocricetus_auratus,Mesoplodon_bidens,Microcebus_murinus,Microgale_talazaci,Micronycteris_hirsuta,Microtus_ochrogaster,Miniopterus_natalensis,Miniopterus_schreibersii,Mirounga_angustirostris,Mirza_coquereli,Monodon_monoceros,Mormoops_blainvillei,Moschus_moschiferus,Mungos_mungo,Murina_feae,Muscardinus_avellanarius,Mus_caroli,Mus_musculus,Mus_pahari,Mus_spretus,Mustela_putorius,Myocastor_coypus,Myotis_brandtii,Myotis_davidii,Myotis_lucifugus,Myotis_myotis,Myrmecophaga_tridactyla,Nannospalax_galili,Nasalis_larvatus,Neomonachus_schauinslandi,Neophocaena_asiaeorientalis,Noctilio_leporinus,Nomascus_leucogenys,Nycticebus_coucang,Ochotona_princeps,Octodon_degus,Odobenus_rosmarus,Odocoileus_virginianus,Okapia_johnstoni,Ondatra_zibethicus,Onychomys_torridus,Orcinus_orca,Orycteropus_afer,Oryctolagus_cuniculus,Otolemur_garnettii,Ovis_aries,Ovis_canadensis,Pan_paniscus,Panthera_onca,Panthera_pardus,Panthera_tigris,Pantholops_hodgsonii,Pan_troglodytes,Papio_anubis,Paradoxurus_hermaphroditus,Perognathus_longimembris,Peromyscus_maniculatus,Petromus_typicus,Phocoena_phocoena,Piliocolobus_tephrosceles,Pipistrellus_pipistrellus,Pithecia_pithecia,Platanista_gangetica,Pongo_abelii,Procavia_capensis,Propithecus_coquereli,Psammomys_obesus,Pteronotus_parnellii,Pteronura_brasiliensis,Pteropus_alecto,Pteropus_vampyrus,Puma_concolor,Pygathrix_nemaeus,Rangifer_tarandus,Rattus_norvegicus,Rhinolophus_sinicus,Rhinopithecus_bieti,Rhinopithecus_roxellana,Rousettus_aegyptiacus,Saguinus_imperator,Saiga_tatarica,Saimiri_boliviensis,Scalopus_aquaticus,Semnopithecus_entellus,Sigmodon_hispidus,Solenodon_paradoxus,Sorex_araneus,Spermophilus_dauricus,Spilogale_gracilis,Suricata_suricatta,Sus_scrofa,Tadarida_brasiliensis,Tamandua_tetradactyla,Tapirus_indicus,Tapirus_terrestris,Thryonomys_swinderianus,Tolypeutes_matacus,Tonatia_saurophila,Tragulus_javanicus,Trichechus_manatus,Tupaia_chinensis,Tupaia_tana,Tursiops_truncatus,Uropsilus_gracilis,Ursus_maritimus,Vicugna_pacos,Vulpes_lagopus,Xerus_inauris,Zalophus_californianus,Zapus_hudsonius,Ziphius_cavirostris,fullTreeAnc239"
options_hal2maf="--refGenome fullTreeAnc239 --onlyOrthologs --targetGenomes $genomes"

# If have to run for chr X and Y: Parameters for X and Y
sexChr_yes=1 # leave on 1 if you have specific sex chr ot calculate 
ancestral_repeat_final_synt_LO=$output_dir"/rmsk."$species_final".Allfam.LO_"$species_reference".bed"
ancestral_repeat_final_AR100kb_hg38X=$output_dir"/rmsk."$species_final".Allfam.LO_"$species_reference"_chrX.100kb.bed" # random 100kb selected on X
maf_100kb_X=$output_dir"/rmsk."$species_final".Allfam.LO_"$species_reference"_chrX.100kb.maf" # maf of random 100kb positions on X
maf_100kbmod_X=$output_dir"/rmsk."$species_final".Allfam.LO_"$species_reference"_chrX.100kb.mod" # model output for chrX
ancestral_repeat_final_AR100kb_hg38Y=$output_dir"/rmsk."$species_final".Allfam.LO_"$species_reference"_chrY.100kb.bed" # random 100kb selected on Y
maf_100kb_Y=$output_dir"/rmsk."$species_final".Allfam.LO_"$species_reference"_chrY.100kb.maf" # maf of random 100kb positions on Y
maf_100kbmod_Y=$output_dir"/rmsk."$species_final".Allfam.LO_"$species_reference"_chrY.100kb.mod" # model output for chrY
options_hal2maf_XY="--refGenome $species_reference --onlyOrthologs --targetGenomes $genomes" # hal2maf parameters

### END OF PARAMETERS

# Extract sequence
hal2fasta $hal_file $species_rmsk > $speciesrmsk_fasta

### RepeatMask the ancestral sequence of reference
RepeatMasker $speciesrmsk_fasta # Output will be in folder with the smae name than ref.

### Filter out problematic repeats # Credits to A.Smit
./exclude.pl $rmskOUT > $rmskOUTF
awk -v OFS='\t' '{print $5,$6,$7}' $rmskOUTF | sort -k1,1 -k2,2n | bedtools merge -i - > $rmskOUTBEDF

### Reconvert to most ancestral sequence : Anc239
if [[ "$species_rmsk" != "$species_final" && "$species_final" != "" ]]; then
	halLiftOver $hal_file $species_rmsk $rmskOUTBEDF $species_final $ancestral_repeat_final --outPSL --noDupes
	sort -k1,1 -k2,2n $ancestral_repeat_final | bedtools merge -i - > $ancestral_repeat_final".temp" && mv $ancestral_repeat_final".temp" $ancestral_repeat_final
fi

### Check synteny with 4 Ancestors of the tree Xenathran, Afrotherian, Laurasitherian and Eurachotonglires
for ancestral_branch in $list_ancestor; do
	outfile=$output_dir"/"$species_final".filtered.merged.vs"$ancestral_branch".fa.bed"
	halLiftover $halfile $species_final $ancestral_repeat_final $ancestral_branch $outfile --noDupes --outPSL
	awk -v OFS='\t' '{print $10,$12,$13,$14,$16,$17}' $outfile > $outfile"_temp" && mv $outfile"_temp" $outfile
	bedtools intersect -a $ancestral_repeat_final_synt -b $outfile > $ancestral_repeat_final_synt"_temp" && mv $ancestral_repeat_final_synt"_temp" $ancestral_repeat_final_synt
done

### Calculate models

# Calculate model on random positions of Anc239 repeat coordinates, alignment referenced to AncRep239
bedops --chop 1 $ancestral_repeat_final_synt | shuf -n 100000 | sort -k1,1 -k2,2n | bedtools merge -i - > $ancestral_repeat_final_AR100kb

# run hal2maf
hal2maf $options_hal2maf --refTargets $ancestral_repeat_final_AR100kb $halfile $maf_100kb
/proj/uppstore2017228/KLT.04.200M/200m_MD/resources/mafTools/bin/mafDuplicateFilter -m $maf_100kb > $maf_100kb"_temp" && mv $maf_100kb"_temp" $maf_100kb

#phyloFit
phyloFit $maf_100kb -i MAF --subst-mod REV --EM --tree $tree --out-root $maf_100kbmod

### Chr X and Y models
# For models on chrX and Y, conversion of coordinates to a species with identified chrX and Y 
# In that case, we picked human.
halLiftOver $hal_file $species_rmsk $ancestral_repeat_final_synt $species_reference $ancestral_repeat_final_synt_LO --outPSL --noDupes

# Select random 100kb on chrX and chrY
awk -v OFS='\t' '{if ($1=="chrX") print}' $ancestral_repeat_final_synt | bedops --chop 1 - | shuf -n 100000 | sort -k1,1 -k2,2n | bedtools merge -i - > $ancestral_repeat_final_AR100kb_hg38X
awk -v OFS='\t' '{if ($1=="chrY") print}' $ancestral_repeat_final_synt | bedops --chop 1 - | shuf -n 100000 | sort -k1,1 -k2,2n | bedtools merge -i - > $ancestral_repeat_final_AR100kb_hg38Y

# Extract / calculate model for X: 
hal2maf $options_hal2maf_XY --refTargets $ancestral_repeat_final_AR100kb_hg38X $hal_file $maf_100kb_X
phyloFit $maf_100kb_X -i MAF --subst-mod REV --EM --tree $tree --out-root $maf_100kbmod_X

# Extract / calculate model for Y:
hal2maf $options_hal2maf_XY --refTargets $ancestral_repeat_final_AR100kb_hg38Y $hal_file $maf_100kb_Y
phyloFit $maf_100kb_Y -i MAF --subst-mod REV --EM --tree $tree --out-root $maf_100kbmod_Y

echo -n "Time ended: "
date
