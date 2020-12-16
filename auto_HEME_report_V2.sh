#!bin/bash

indir=/data/IR/data/IR_Org/ion.reporter@lifetech.com
outdir=/mnt/Z_drive/acc_pathology/molecular/MOLECULAR/IonTorrent/HEME_Report/$1

mkdir -p ${outdir}
mkdir -p ${outdir}/API_zip
mkdir -p ${outdir}/Report

Rscript samplelist_create_V2.R $1 $2

###### to call API (token generated on 09/14/20  https://10.182.97.92/ir/secure/home.html)
curl -v -k -X GET "https://10.182.97.92/api/v1/getvcf?format=json&start_date=$3&end_date=$4" -H "Content-Type:application/x-www-form-urlencoded" -H "Authorization:YTYzMTI5YWY0MTY3Y2FhNzdkMmNkOTgxMWNmZjJiODc0NGNiN2U3MWVmMjg4MGUwZTFjY2YwYWJjZDc2MjBjYQ"

echo " "
echo API call succeeded on Status Code 200, otherwise need to re-call API.

for i in $(cat ${outdir}/${1}_samplelist_V2.txt | awk '{print $1}')
do
	sam_dir=$(ls -lt $indir | grep $i | head -n1)  
	#echo $sam_dir
	sam_dir_single=$(echo $sam_dir |awk '{print $9}')
	#echo $sam_dir_single
	if [ -z $sam_dir_single ]
	then
		echo $i on the wetlab spreedsheet is NOT FOUND on IR backend. Need manually check wetlab spreedsheet 'in' /MOLECULAR LAB ONLY/Oncomine Patient Data/Worksheets/${1}.xlse
		echo " "
	elif [ ! -z $sam_dir_single ]
	then
		echo $sam_dir_single >> ${outdir}/${1}_sampledir_IR.txt
	echo $sam_dir_single is written in new text file : ${outdir}/${1}_sampledir_IR.txt for further processing data     	
	fi
done

cat ${outdir}/${1}_sampledir_IR.txt | sort |uniq -d >> ${outdir}/${1}_sampledir_IR_final.txt
cat ${outdir}/${1}_sampledir_IR.txt | sort |uniq -u >> ${outdir}/${1}_sampledir_IR_final.txt

for j in $(cat ${outdir}/${1}_sampledir_IR_final.txt | awk '{print $1}')
do
	sam_id=$(echo $j|cut -d'_' -f 1)
	echo $sam_id
	zip_dir=$(ls -lt ${indir}/${j} | grep $sam_id | head -n1| awk '{print $9}')
	zip_file=$(ls -lt ${indir}/${j}/${zip_dir}/*.zip | head -n1 | awk '{print $9}')
	sudo scp ${zip_file} ${outdir}/API_zip
	echo ${zip_file_single} has been copied to ${outdir}/API_zip
	echo " "
done

for z in $(ls ${outdir}/API_zip/*.zip)
do
    sample=$(echo $z | cut -d'/' -f 7 | cut -d'_' -f 1-2)
    unzip -d ${outdir} $z
    echo $z
done

Rscript HEME_full_report_V3 $1 $2
