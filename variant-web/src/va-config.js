CELLBASE_HOST = "http://ws.bioinfo.cipf.es/cellbase/rest";
CELLBASE_VERSION = "v3";
OPENCGA_HOST = "http://ws-beta.bioinfo.cipf.es/opencga/rest";

if (
    window.location.host.indexOf("localhost") != -1 ||
        window.location.host.indexOf("rsanchez") != -1 ||
        window.location.host.indexOf("imedina") != -1 ||
        window.location.href.indexOf("http://bioinfo.cipf.es/apps-beta") != -1 ||
        window.location.protocol === "file:"
    ) {
    CELLBASE_HOST = "http://ws-beta.bioinfo.cipf.es/cellbase/rest";
    OPENCGA_HOST = "http://ws-beta.bioinfo.cipf.es/opencga-staging/rest";
}

var h = window.location.host;
if (h.indexOf("fsalavert") != -1) {
    CELLBASE_HOST = "http://ws-beta.bioinfo.cipf.es/cellbase/rest";
    OPENCGA_HOST = "http://ws-beta.bioinfo.cipf.es/opencga-staging-paco/rest";
}
if (h.indexOf("aaleman") != -1) {
    CELLBASE_HOST = "http://ws-beta.bioinfo.cipf.es/cellbase/rest";
    OPENCGA_HOST = "http://ws-beta.bioinfo.cipf.es/opencga-staging-aleman/rest";
}


SUITE_INFO = '<div style=" width: 800px;">'
    + '<h1 style="font-family: Source Sans Pro;">Overview</h1><span align="justify">VARIANT (VARIant ANalysis Tool) can report the functional properties of any variant in all the human, mouse or rat genes (and soon new model organisms will be added) and the corresponding neighborhoods. Also other non-coding extra-genic regions, such as miRNAs are included in the analysis.<br><br>	VARIANT not only reports the obvious functional effects in the coding regions but also analyzes noncoding SNVs situated both within the gene and in the neighborhood that could affect different regulatory motifs, splicing signals, and other structural elements. These include: Jaspar regulatory motifs, miRNA targets, splice sites, exonic splicing silencers, calculations of selective pressures on the particular polymorphic positions, etc.</span>'
    + '<br><br><br>'
    + '<p align="justify"><h2 style="font-family: Source Sans Pro;">Note</h2>This web application makes an intensive use of new web technologies and standards like HTML5, so browsers that are fully supported for this site are: Chrome 14+, Firefox 7+, Safari 5+ and Opera 11+. Older browser like Chrome13-, Firefox 5- or Internet Explorer 9 may rise some errors. Internet Explorer 6 and 7 are no supported at all.</p>'
    + '</div>';
//    +'</div>'+
//    +'<br><br><h2>Sign in</h2><p style=" width: 800px;">You must be logged in to use this Web application, you can <b><i>register</i></b> or use a <b><i>anonymous user</i></b> as shown in the following image by clicking on the <b><i>"Sign in"</i></b> button on the top bar</p><br><div style="float:left;"><img src="http://jsapi.bioinfo.cipf.es/libs/resources/img/loginhelpbutton.png"></div><img src="http://jsapi.bioinfo.cipf.es/libs/resources/img/loginhelp.png">';
//    +'';


var POPULAR_SPECIES = ["Homo sapiens", "Mus musculus", "Danio rerio", "Drosophila melanogaster", "Saccharomyces cerevisiae", "Plasmodium falciparum", "Arabidopsis thaliana", "Citrus Clementina"];

var AVAILABLE_SPECIES = {
    "text": "Species",
    "items": [
        {
            "text": "Vertebrates",
            "items": [
                {"text": "Bos taurus", "assembly": "UMD3.1", "region": {"chromosome": "1", "start": 1000000, "end": 1000000}, "chromosomes": ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "X", "MT"], "url": "ftp://ftp.ensembl.org/pub/release-71/"},
                {"text": "Canis familiaris", "assembly": "CanFam3.1", "region": {"chromosome": "1", "start": 1000000, "end": 1000000}, "chromosomes": ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "X", "MT"], "url": "ftp://ftp.ensembl.org/pub/release-71/"},
                {"text": "Ciona intestinalis", "assembly": "KH", "region": {"chromosome": "1", "start": 1000000, "end": 1000000}, "chromosomes": ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "MT"], "url": "ftp://ftp.ensembl.org/pub/release-71/"},
                {"text": "Danio rerio", "assembly": "Zv9", "region": {"chromosome": "1", "start": 1000000, "end": 1000000}, "chromosomes": ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "X", "Y", "MT"], "url": "ftp://ftp.ensembl.org/pub/release-71/"},
                {"text": "Equus caballus", "assembly": "EquCab2", "region": {"chromosome": "1", "start": 1000000, "end": 1000000}, "chromosomes": ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "X", "MT"], "url": "ftp://ftp.ensembl.org/pub/release-71/"},
                {"text": "Felis catus", "assembly": "Felis_catus_6.2", "region": {"chromosome": "1", "start": 1000000, "end": 1000000}, "chromosomes": ["A1", "A2", "A3", "B1", "B2", "B3", "B4", "C1", "C2", "D1", "D2", "D3", "D4", "E1", "E2", "E3", "F1", "F2", "X", "MT"], "url": "ftp://ftp.ensembl.org/pub/release-71/"},
                {"text": "Gallus gallus", "assembly": "Galgal4", "region": {"chromosome": "1", "start": 1000000, "end": 1000000}, "chromosomes": ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "32", "W", "Z", "MT"], "url": "ftp://ftp.ensembl.org/pub/release-71/"},
                {"text": "Gorilla gorilla", "assembly": "gorGor3.1", "region": {"chromosome": "1", "start": 1000000, "end": 1000000}, "chromosomes": ["1", "2a", "2b", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "MT"], "url": "ftp://ftp.ensembl.org/pub/release-71/"},
                {"text": "Homo sapiens", "assembly": "GRCh37.p10", "region": {"chromosome": "13", "start": 32889611, "end": 32973805}, "chromosomes": ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT"], "url": "ftp://ftp.ensembl.org/pub/release-71/"},
//                {"text": "Homo sapiens", "assembly": "GRCh37.p10", "region": {"chromosome": "13", "start": 32889599, "end": 32889739}, "chromosomes": ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT"], "url": "ftp://ftp.ensembl.org/pub/release-71/"},
//                {"text": "Homo sapiens", "assembly": "GRCh37.p10", "region": {"chromosome": "20", "start": 32878277, "end": 32878277}, "chromosomes": ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT"], "url": "ftp://ftp.ensembl.org/pub/release-71/"},
//                {"text": "Homo sapiens", "assembly": "GRCh37.p10", "region": {"chromosome": "1", "start": 32877109, "end": 32882337}, "chromosomes": ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT"], "url": "ftp://ftp.ensembl.org/pub/release-71/"},
                {"text": "Macaca mulatta", "assembly": "MMUL_1.0", "region": {"chromosome": "1", "start": 1000000, "end": 1000000}, "chromosomes": ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "X", "MT"], "url": "ftp://ftp.ensembl.org/pub/release-71/"},
                {"text": "Mus musculus", "assembly": "GRCm38.p1", "region": {"chromosome": "1", "start": 18422009, "end": 18422009}, "chromosomes": ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "X", "Y", "MT"], "url": "ftp://ftp.ensembl.org/pub/release-71/"},
                {"text": "Oryctolagus cuniculus", "assembly": "oryCun2", "region": {"chromosome": "1", "start": 1000000, "end": 1000000}, "chromosomes": ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "X", "MT"], "url": "ftp://ftp.ensembl.org/pub/release-71/"},
                {"text": "Pan troglodytes", "assembly": "CHIMP2.1.4", "region": {"chromosome": "1", "start": 1000000, "end": 1000000}, "chromosomes": ["1", "2A", "2B", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT"], "url": "ftp://ftp.ensembl.org/pub/release-71/"},
                {"text": "Pongo abelii", "assembly": "PPYG2", "region": {"chromosome": "1", "start": 1000000, "end": 1000000}, "chromosomes": ["1", "2a", "2b", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Un", "MT"], "url": "ftp://ftp.ensembl.org/pub/release-71/"},
                {"text": "Rattus norvegicus", "assembly": "Rnor_5.0", "region": {"chromosome": "1", "start": 1000000, "end": 1000000}, "chromosomes": ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "X", "MT"], "url": "ftp://ftp.ensembl.org/pub/release-71/"},
                {"text": "Sus scrofa", "assembly": "Sscrofa10.2", "region": {"chromosome": "1", "start": 1000000, "end": 1000000}, "chromosomes": ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "X", "Y", "MT"], "url": "ftp://ftp.ensembl.org/pub/release-71/"}
            ]
        },
        {
            "text": "Metazoa",
            "items": [
                {"text": "Anopheles gambiae", "assembly": "AgamP3", "region": {"chromosome": "2L", "start": 1000000, "end": 1000000}, "chromosomes": ["2L", "2R", "3L", "3R", "X"], "url": "ftp://ftp.ensemblgenomes.org/pub/metazoa/release-18/"},
                {"text": "Caenorhabditis elegans", "assembly": "WBcel235", "region": {"chromosome": "I", "start": 1000000, "end": 1000000}, "chromosomes": ["I", "II", "III", "IV", "V", "X", "MtDNA"], "url": "ftp://ftp.ensemblgenomes.org/pub/metazoa/release-18/fasta/"},
                {"text": "Drosophila melanogaster", "assembly": "BDGP5", "region": {"chromosome": "2L", "start": 1000000, "end": 1000000}, "chromosomes": ["2L", "2LHet", "2R", "2RHet", "3L", "3LHet", "3R", "3RHet", "4", "U", "Uextra", "X", "XHet", "YHet", "dmel_mitochondrion_genome"], "url": "ftp://ftp.ensemblgenomes.org/pub/metazoa/release-18/"},
                {"text": "Drosophila simulans", "assembly": "WUGSC1", "region": {"chromosome": "2L", "start": 1000000, "end": 1000000}, "chromosomes": ["2L", "2R", "3L", "3R", "4", "X"], "url": "ftp://ftp.ensemblgenomes.org/pub/metazoa/release-18/"},
                {"text": "Drosophila yakuba", "assembly": "dyak_r1.3", "region": {"chromosome": "2L", "start": 1000000, "end": 1000000}, "chromosomes": ["2L", "2R", "3L", "3R", "4", "chr2h", "chr3h", "chrXh", "chrYh", "X"], "url": "ftp://ftp.ensemblgenomes.org/pub/metazoa/release-18/"}
            ]
        },
        {
            "text": "Fungi",
            "items": [
                {"text": "Aspergillus fumigatus", "assembly": "CADRE", "region": {"chromosome": "I", "start": 1000000, "end": 1000000}, "chromosomes": ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "MT"], "url": "ftp://ftp.ensemblgenomes.org/pub/fungi/release-18/"},
                {"text": "Aspergillus nidulans", "assembly": "ASM1142v1", "region": {"chromosome": "I", "start": 1000000, "end": 1000000}, "chromosomes": ["I", "II", "III", "IV", "V", "VI", "VII", "VIII"], "url": "ftp://ftp.ensemblgenomes.org/pub/fungi/release-18/"},
                {"text": "Aspergillus niger", "assembly": "CADRE", "region": {"chromosome": "I", "start": 1000000, "end": 1000000}, "chromosomes": ["I", "II", "III", "IV", "V", "VI", "VII", "VIII"], "url": "ftp://ftp.ensemblgenomes.org/pub/fungi/release-18/"},
                {"text": "Aspergillus oryzae", "assembly": "CADRE2", "region": {"chromosome": "I", "start": 1000000, "end": 1000000}, "chromosomes": ["I", "II", "III", "IV", "V", "VI", "VII", "VIII"], "url": "ftp://ftp.ensemblgenomes.org/pub/fungi/release-18/"},
                {"text": "Saccharomyces cerevisiae", "assembly": "SacCer_Apr2011", "region": {"chromosome": "I", "start": 1000000, "end": 1000000}, "chromosomes": ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI", "Mito"], "url": "ftp://ftp.ensemblgenomes.org/pub/fungi/release-18/"},
                {"text": "Schizosaccharomyces pombe", "assembly": "ASM294v1", "region": {"chromosome": "I", "start": 1000000, "end": 1000000}, "chromosomes": ["AB325691", "I", "II", "III", "MT", "MTR"], "url": "ftp://ftp.ensemblgenomes.org/pub/fungi/release-18/"}
            ]
        },
        {
            "text": "Protist",
            "items": [
                {"text": "Leishmania major", "assembly": "ASM272v2", "region": {"chromosome": "1", "start": 1000000, "end": 1000000}, "chromosomes": ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36"], "url": "ftp://ftp.ensemblgenomes.org/pub/protists/release-18/"},
                {"text": "Plasmodium falciparum", "assembly": "ASM276v1", "region": {"chromosome": "01", "start": 1000000, "end": 1000000}, "chromosomes": ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14"], "url": "ftp://ftp.ensemblgenomes.org/pub/protists/release-18/"}
            ]
        },
        {
            "text": "Plants",
            "items": [
                {"text": "Arabidopsis lyrata", "assembly": "v.1.0", "region": {"chromosome": "1", "start": 1000000, "end": 1000000}, "chromosomes": ["1", "2", "3", "4", "5", "6", "7", "8"], "url": "ftp://ftp.ensemblgenomes.org/pub/plants/release-18/"},
                {"text": "Arabidopsis thaliana", "assembly": "TAIR10", "region": {"chromosome": "1", "start": 1000000, "end": 1000000}, "chromosomes": ["1", "2", "3", "4", "5", "Mt", "Pt"], "url": "ftp://ftp.ensemblgenomes.org/pub/plants/release-18/"},
                {"text": "Brachypodium distachyon", "assembly": "v1.0", "region": {"chromosome": "1", "start": 1000000, "end": 1000000}, "chromosomes": ["1", "2", "3", "4", "5"], "url": "ftp://ftp.ensemblgenomes.org/pub/plants/release-18/"},
                {"text": "Glycine max", "assembly": "V1.0", "region": {"chromosome": "Gm01", "start": 1000000, "end": 1000000}, "chromosomes": ["Gm01", "Gm02", "Gm03", "Gm04", "Gm05", "Gm06", "Gm07", "Gm08", "Gm09", "Gm10", "Gm11", "Gm12", "Gm13", "Gm14", "Gm15", "Gm16", "Gm17", "Gm18", "Gm19", "Gm20"], "url": "ftp://ftp.ensemblgenomes.org/pub/plants/release-18/"},
                {"text": "Oryza sativa", "assembly": "MSU6", "region": {"chromosome": "1", "start": 1000000, "end": 1000000}, "chromosomes": ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "Mt", "Pt", "Sy", "Un"], "url": "ftp://ftp.ensemblgenomes.org/pub/plants/release-18/"},
                {"text": "Vitis vinifera", "assembly": "IGGP_12x", "region": {"chromosome": "1", "start": 1000000, "end": 1000000}, "chromosomes": ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "Un"], "url": "ftp://ftp.ensemblgenomes.org/pub/plants/release-18/"},
                {"text": "Zea mays", "assembly": "AGPv3", "region": {"chromosome": "1", "start": 1000000, "end": 1000000}, "chromosomes": ["1", "2", "3", "4", "5", "6", "7", "8", "9", "Mt", "Pt"], "url": "ftp://ftp.ensemblgenomes.org/pub/plants/release-18/"}
//                  {"text": "Citrus Clementina", "assembly": "1.0", "region":{"chromosome":"scaffold_1","start":233423,"end":236969},   "chromosomes": [], "url": ""}
            ]
        }
    ]
};


/** Reference to a species from the list to be shown at start **/
var DEFAULT_SPECIES = AVAILABLE_SPECIES.items[0].items[8];
