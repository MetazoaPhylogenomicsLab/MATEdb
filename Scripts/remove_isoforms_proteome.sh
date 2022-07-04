#!usr/bin/bash

############################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo "Script created by Gemma I. Martinez-Redondo to obtain the longest isoform of each gene from a protein FASTA file of a genome."
   echo
   echo "Syntax: bash remove_isoforms_proteome.sh [-g|f|s|o|h]"
   echo "options:"
   echo "-g/--gff           Input GFF file. Take into account that the content of the last GFF column is not standarized. Check GFF file and line 128 of this script to confirm that the script will work with your GFF file. Regex are used."
   echo "-h/--help             Print this Help."
   echo "-f/--fasta          Input protein FASTA file. Sequences must be sorted alphabetically by header."
   echo "-s/--species             4-letter+1-number code of the species."
   echo "-o/--filtered-output     (Optional) Output file name. If no ouptut file name is given, 'inputfile.filtered.fasta' will be used."
   echo "-t/--trinity-header-proteome (Optional) Path to trinity_header_proteome.py script. If not provided, current directory will be used as default."
   echo "-l/--fetch-longest (Optional) Path to fetch_longest_iso.py script. If not provided, current directory will be used as default."
   echo
}

############################################################
############################################################
# Extra functions (print error message)                    #
############################################################
############################################################
die() { echo "$*" 1>&2 ; exit 1; }

############################################################
############################################################
# Check arguments                                          #
############################################################
############################################################

#Check argument options.
while :; do
    case $1 in
        -g|--gff)
            #Check if a GFF file has been given.
        if [[ "$2" ]]; then
            GFF=$2
            shift
        else
            die 'ERROR: "--gff" requires a non-empty option argument.'
        fi
        ;;
        -f|--fasta)
            #Check if a protein fasta file has been given.
            if [[ "$2" ]]; then
            PROTS=$2
            shift
        else
            die 'ERROR: "--fasta" requires a non-empty option argument.'
        fi
        ;;
        -s|--species)
            #Check if a species code has been given.
            if [[ "$2" ]]; then
            SPECIES=$2
            shift
        else
            die 'ERROR: "--species" requires a non-empty option argument.'
        fi
        ;;
        -o|--filtered-output)
            #Check if an output file has been given.
            if [[ "$2" ]]; then
            OUTFILE=$2
            shift
        else
            die 'ERROR: "--filtered-output" requires a non-empty option argument.'
        fi
        ;;
	-t|--trinity-header-proteome)
	    #Check if path to trinity_header_proteome.py
	    if [[ "$2" ]]; then
	    	TRINITY_HEADER_PROTEOME=$2
	    shift
	;;
	-l|--fetch-longest)
            #Check if path to fetch_longest_iso.py
            if [[ "$2" ]]; then
                FETCH_LONGEST=$2
            shift
        ;;
        -h|--help)
            #Print help.
            Help
            exit
        ;;
        *) break
    esac
    shift
done

#If no output file name have been given, use GFF file name.
if [[ ! "$OUTFILE" ]]; then
    OUTPATH=$(basename ${GENETREE%.*})
    OUTFILE=${GFF%.*}.filtered.fasta
fi

#If no fetch_longest_iso.py or trinity_header_proteome.py path have been given, use current directory
if [[ ! "$TRINITY_HEADER_PROTEOME" ]]; then
    TRINITY_HEADER_PROTEOME=trinity_header_proteome.py
fi
if [[ ! "$FETCH_LONGEST" ]]; then
    FETCH_LONGEST=fetch_longest_iso.py
fi

############################################################
############################################################
# Main program                                             #
############################################################
############################################################

#Remove isoforms from downloaded proteomes.
#Obtain CDS from GFF file and fuse the ones that have the same gene (second column of geneprot.txt file) and protein id (first column). These CDS correspond to the same isoform and are redundant.
#WARNING!
#THIS STEP REQUIRES KNOWING THE CONTENT OF THE LAST GFF COLUMN. STRUCTURE USUALLY FOUND IN NCBI GENOMES IS WRITTEN, BUT THIS CAN VARY
grep -vE "^#" $GFF \
| awk 'BEGIN{FS="\t"} {if($3=="CDS") {print $9}}' \
| sed 's/ID=cds-\([^\;]*\)\;*.*gene=\([^\;]*\).*/\1 \2/g' \
| sort | uniq > ${SPECIES}_geneprot.txt

#Check that number of proteins is correct.
if [[ $(wc -l ${SPECIES}_geneprot.txt | cut -d" " -f1) != $(grep -c ">" $PROTS) ]]
then
	echo "Number of proteins differs between annotation and proteome file"
	exit
fi

#Important! Both files must have genes in the same order.
cut -d" " -f1 ${SPECIES}_geneprot.txt > order_${SPECIES}_geneprot.txt
grep ">" $PROTS | cut -d" " -f1 | cut -d">" -f2 > order_${SPECIES}.fasta
if [[ $(diff order_${SPECIES}_geneprot.txt order_${SPECIES}.fasta) != "" ]]
then
	echo "Proteins in the proteome file are not sorted"
	exit
fi
rm {order_${SPECIES}_geneprot.txt,order_${SPECIES}.fasta}

#Save original headers.
grep ">" $PROTS >${SPECIES}_orig_headers.txt

#Execute python script for changing headers to Trinity-style headers
python $TRINITY_HEADER_PROTEOME -i ${SPECIES}_geneprot.txt -o ${SPECIES}_headers.txt -s $SPECIES

#Change headers with the new ones
< ${SPECIES}_headers.txt  perl -pe '$_ = <STDIN> if /^>/' $PROTS >${SPECIES}.mod.fasta

#Save new and old headers into a conversion file for further use
paste ${SPECIES}_orig_headers.txt ${SPECIES}_headers.txt > ${SPECIES}_conversion.txt
rm {${SPECIES}_headers.txt,${SPECIES}_orig_headers.txt}

#Remove isoforms as usual
python $FETCH_LONGEST -i ${SPECIES}.mod.fasta -o ${SPECIES}.filtered.fasta -t -l
