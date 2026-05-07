Working stationary test:

./gps-sim -e brdc0980.26n -l=44.2363968,-71.0091981,10 -d 300 -r iqfile --iq16 -I

sudo sysctl -w net.core.rmem_max=2500000

sudo sysctl -w net.core.wmem_max=2500000

/usr/local/lib/uhd/examples/tx_samples_from_file --file iqdata.bin --type short --rate 3e6 --freq 1575.42e6 --gain 100 --repeat --args "addr=192.168.1.230"

Working motion file circle.csv:

./gps-sim -e brdc0980.26n -m circle.csv -d 300 -r iqfile --iq16 -I

/usr/local/lib/uhd/examples/tx_samples_from_file --file iqdata.bin --type short --rate 3e6 --freq 1575.42e6 --gain 100 --repeat --args "addr=192.168.1.230"

LAT: 35.27N

LON: 137.03E

gpspipe -rt | tee 12500M-test.out | grep GPGGA

cgps --silent
