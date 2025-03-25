# downloading the data from CSC bucket via ssh client on LINUX
# https://docs.csc.fi/data/Allas/using_allas/rclone_local/

sudo -v ; curl https://rclone.org/install.sh | sudo bash # install rclone
wget https://raw.githubusercontent.com/CSCfi/allas-cli-utils/master/allas_conf
source allas_conf -u CSCUSERNAME -p PROJECTNAME
rclone lsd allas: # list buckets within project
rclone lsd allas:BUCKET # list files within bucket
rclone copy allas:BUCKET LOCAL/PATH # copy files from csc bucket to the local directory - may be done the other way as well
rclone sync --interactive allas:BUCKET /LOCAL/PATH  # syncs csc bucket to the local - may be done the other way as well


# encryption-decryption manual
# https://docs.csc.fi/data/sensitive-data/sequencing_center_tutorial/


