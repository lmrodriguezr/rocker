# ROCker

Accurately detecting functional genes in metagenomes.

## System requirements

1. [Ruby](https://www.ruby-lang.org/), with the [restclient](https://rubygems.org/gems/rest_client) and
   [nokogiri](http://www.nokogiri.org/) packages. To install the required packages execute:
   
   ```bash
   $> gem install rest_client
   $> gem install nokogiri
   $> gem install json
   ```

2. [R](http://www.r-project.org/), with the [pROC](http://cran.r-project.org/web/packages/pROC/index.html)
   package. To install the required packages, execute:

   ```bash
   $> R
   R> install.packages('pROC');
   ```

3. A metagenome simulation software: [Grinder](http://sourceforge.net/projects/biogrinder/).

4. A local search software: [NCBI BLAST+](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/),
   [DIAMOND](http://ab.inf.uni-tuebingen.de/software/diamond/), or any other software producing a results
   in the same format of Tabular-BLAST without comments (use the `--search-cmd` and `--makedb-cmd` options
   to execute any other software).

5. A multiple alignment software: [MUSCLE](http://www.drive5.com/muscle/),
   [Clustal Omega](http://www.clustal.org/omega/) or any other software supporting FastA input and output
   (use the `--aligner-cmd` option to execute any other software).

## Installation

Install ROCker using [RubyGems](https://rubygems.org/gems/bio-rocker):

```bash
$> gem install bio-rocker
```

Or get the source from [GitHub](https://github.com/lmrodriguezr/rocker):

```bash
$> git clone https://github.com/lmrodriguezr/rocker.git
$> ./rocker/bin/rocker
```

## License

[Artistic license 2.0](http://www.perlfoundation.org/artistic_license_2_0).

## Authors

Luis H (Coto) Orellana, Luis M. Rodriguez-R & Konstantinos Konstantinidis, at the
[Kostas lab](http://enve-omics.gatech.edu/).

