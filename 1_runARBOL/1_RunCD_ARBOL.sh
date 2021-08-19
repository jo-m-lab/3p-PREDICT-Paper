export $(cat .env | xargs)
dsub \
  --provider google-v2 \
  --project shalek-lab \
  --logging $CD_BUCKET_FOLDER/logs \
  --regions us-central1 \
  --boot-disk-size 15  --timeout 7d \
  --min-ram 256 --min-cores 4 --disk-size 500 \
  --input STARTING_SROBJ=$CD_BUCKET_FOLDER/CDmerged_srobj.rds \
  --output-recursive OUT=$CD_BUCKET_FOLDER/$(date +%Y%m%d)_ITCoutput \
  --image shaleklab/terra-seurat-env:latest \
  --script GenTeiredClusters.R