Reporting pipeline for cancer diagnostic NGS panel on ThermoFisher IonTorrent platform

Steps to run script on IR 
1.	Log in IR (Ion Report) backend
2.	cd /home/ionadmin/script (where script saved)
3.	sh auto_HEME_report_V2.sh $1 $2 $3 $4| tee heme_run_id.log (e.g.heme_20-xx.log)
4.	output report will be saved in $2/$1/Report/run_id_heme_report.csv

Notes: 
$1: run_id (e.g.HEME-20-xx)
$2: target_dir e.g /mnt/Z_drive/acc_pathology/molecular/MOLECULAR/IonTorrent/HEME_report/
$3: sequencing start date (yyyy-mm-dd)
$4: sequencing end date (yyyy-mm-dd)

Example: 
sh auto_HEME_report_V2.sh heme_20-25 /mnt/Z_drive/acc_pathology/molecular/MOLECULAR/IonTorrent/HEME_report 2020-11-01 2020-11-02