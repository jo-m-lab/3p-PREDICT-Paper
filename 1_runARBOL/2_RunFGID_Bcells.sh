export $(cat .env | xargs)
dsub \
  --provider google-v2 \
  --project shalek-lab \
  --logging $FG_BUCKET_FOLDER/logs \
  --regions us-central1 \
  --boot-disk-size 15  --timeout 7d \
  --min-ram 256 --min-cores 4 --disk-size 500 \
  --input STARTING_SROBJ=$FG_BUCKET_FOLDER/output/srobjs/tier1_cluster1/subset_srobj.rds \
  --output-recursive OUT=$FG_BUCKET_FOLDER/$(date +%Y%m%d)_ITC_BCells_output \
  --image shaleklab/terra-seurat-env:latest \
  --script GenTeiredClustersBcells.R



