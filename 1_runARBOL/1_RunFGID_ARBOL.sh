export $(cat .env | xargs)
dsub \
  --provider google-v2 \
  --project shalek-lab \
  --logging $FG_BUCKET_FOLDER/logs \
  --regions us-central1 \
  --boot-disk-size 15  --timeout 7d \
  --min-ram 256 --min-cores 4 --disk-size 500 \
  --input STARTING_SROBJ=$FG_BUCKET_FOLDER/FGfull_srobj.rds \
  --output-recursive OUT=$FG_BUCKET_FOLDER/output \
  --image shaleklab/terra-seurat-env:latest \
  --script GenTeiredClusters.R