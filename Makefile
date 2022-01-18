libhts.a: Dependencies/htslib-1.9
	cd Dependencies/htslib-1.9; chmod +x ./configure ./version.sh
	cd Dependencies/htslib-1.9; ./configure CFLAGS=-fpic --disable-bz2 --disable-lzma  && make

gpu: requirements.txt libhts.a generate.cpp models.cpp gen.cpp setup.py
	pip install -r requirements.txt;
	python3 setup.py build_ext; python3 setup.py install
	rm -rf build

cpu: requirements_cpu.txt libhts.a generate.cpp models.cpp gen.cpp setup.py
	pip install -r requirements_cpu.txt; pip install torch==1.3.1+cpu torchvision==0.4.2+cpu -f https://download.pytorch.org/whl/torch_stable.html
	python3 setup.py build_ext; python3 setup.py install
	rm -rf build

clean:
	cd Dependencies/htslib-1.9 && make clean || exit 0


