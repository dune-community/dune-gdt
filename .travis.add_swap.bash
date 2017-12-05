#/bin/bash
set -e
SWAP=$HOME/swap.img
MB_SIZE=${1:-4000}
sudo -E dd if=/dev/zero of=${SWAP} bs=1M count=${MB_SIZE}
sudo -E chown root:root ${SWAP}
sudo -E chmod 0600 ${SWAP}
sudo -E mkswap ${SWAP}
sudo -E swapon ${SWAP}
echo "echo 3 > /proc/sys/vm/drop_caches" | sudo -E sh

