class Settings(object):
    def __init__(self):
        self.format = "fasta"
        self.species = 'human'
        self.inputfile = None
        self.outputDir= None
        self.overlapLen =  5
        self.overlapStep = False
        self.filterThreshold = 1
        self.extendedOutput = False
        self.is_digGold = False
        self.noCast = False
        self.castThreshold = {'IGH': 0.2,'IGK': 0.2,'IGL': 0.2, 'TRA': 0.3, 'TRB': 0.2, 'TRD': 0.2, 'TRG': 0.2}
        self.chains = ['IGH', 'IGK', 'IGL', 'TRA', 'TRB', 'TRD', 'TRG']
        self.minlen1 = 2
        self.minlen2 = 1
        self.mismatch1 = 2
        self.mismatch2 = 2
        self.sampleName = None
        self.kmer_len = 3
        self.allReads = None
        self.is_hg38 = False
        self.is_chrFormat2 = False



    def configure(self, **kwargs):
        for name, value in kwargs.items():
            setattr(self, name, value)


    def checkValues(self):
        #Check compatibility of options
        if args.is_hg38 and args.species=="mouse":
            print("::::::ERROR. Options --hg38 and -s mouse are not compatible. Please keep only one of those options.")
            print("Exit!")
            sys.exit(1)

        if args.format == "bam" and args.is_digGold:
            print("::::::ERROR. Options --bam and --digGold are not compatible. Please keep only one of those options.")
            sys.exit(1)

        if args.species:
            if args.species in ["human", "mouse"]:
                set_dict["species"] = args.species
            else:
                raise Exception("Species must be either human or mouse")
