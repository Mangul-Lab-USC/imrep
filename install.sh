pip install pysam --user
pip install biopython==1.76 --user
pip install intervaltree --user
pip install https://files.pythonhosted.org/packages/cb/bb/2362099ca5d680e39f75a37b2c8f677fc2d3dda94ce51b3738feff58d136/jellyfish-0.6.0.tar.gz#sha256=f5da646c4ff578ff37695915c05eaae938b1218082bfef59f8a7183477681ab0 --user
pip install numpy --user
pip install networkx --user
git clone https://github.com/jaquejbrito/suffix_tree.git
cd suffix_tree
python setup.py install --user
cd ..
rm -rf suffix_tree
