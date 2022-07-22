# install ocean tracker
# run from inside dir where install required

rm oceantracker -rf # replacing any existing package
git clone https://github.com/oceantracker/oceantracker.git
cd ./oceantracker
python3 -m venv venv
source ./venv/bin/activate
python setup.py develop
pip install -r ./requirements.txt

# test by running one demo
echo "Test run demo 60"
cd ./demos
python run_demos.py -n --d 2
cd ..
cd ..







