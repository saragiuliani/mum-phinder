
wget https://ftp.ebi.ac.uk/pub/databases/covid19dataportal/cdp-file-downloader.zip
unzip cdp-file-downloader.zip

mkdir -p ${1}

for i in {00..39};
do
    acc="accessions/acc_${i}"
    java -jar cdp-file-downloader.jar --domain=VIRAL_SEQUENCES --datatype=SEQUENCES --format=FASTA --location=${1} --email=NONE --accessions=`awk '{gsub("\"",""); print $1}' ${acc} | paste -s -d, -` --protocol=FTP --asperaLocation=null
done

touch ${1}/SARS-CoV2.400k.fa
files=(`ls ${1}/viral_sequences/sequences/fasta/`)
for i in ${files[@]};
do
    cat ${1}/viral_sequences/sequences/fasta/${i} >> ${1}/SARS-CoV2.400k.fa
    rm ${1}/viral_sequences/sequences/fasta/${i}
done

wget -O ${1}/MZ477765.fa https://www.ebi.ac.uk/ena/browser/api/fasta/MZ477765.1?download=true