hs:
	$(MAKE) -C hs lib install

clean:
	rm -vf *~
	$(foreach i, $(shell ls -d */), $(MAKE) -C $(i) clean;)

.PHONY: hs clean
