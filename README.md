# ROCker

Accurately detecting functional genes in metagenomes.

For installation instructions, see [INSTALL.md](./INSTALL.md). For license information, see
[LICENSE.txt](./LICENSE.txt).

## Using existing models

Once you have installed ROCker, the easiest way to use it is by searching pre-existing
models. We maintain a [list of precomputed models](http://enve-omics.ce.gatech.edu/rocker/models)
that you're free to use.

1. **Obtain the model** of interest either downloading it from our repository or creating one yourself
   (see below).

2. **Execute ROCker search**. The minimum required parameters are:

   ```bash
   $> ROCker search -q input.fasta -k model.rocker -o output.blast 
   ```

   Where `input.fasta` is the input metagenome in FastA format, `model.rocker` is the ROCker model,
   and `output.blast` is the output file to be created in tabular BLAST format. For additional
   supported options, execute `ROCker search -h`.

3. **If you have a pre-computed BLAST file**, you can execute instead:

   ```bash
   $> ROCker filter -x input.blast -k model.rocker -o output.blast 
   ```
   
   Where `input.blast` is the input search to be filtered in tabular BLAST format, `model.rocker` is the
   ROCker model, and `output.blast` is the output file to be created in tabular BLAST format. For additional
   supported options, execute `ROCker filter -h`.

## Creating models

Collect a good reference collection of the gene of interest. This is the most important step, but there are
some resources to help you. In general, we find the resources at [UniProt](http://uniprot.org/) very useful.

1. Create a list of **UniProt identifiers** (IDs and/or accessions) representing proteins of the family of
   interest, in a raw text file (one per line).

2. If you want to explicitly exclude certain proteins from the model (*e.g.*, if there are very similar proteins
   with distinct functional properties), create a similar list with those, we will refer to them as a **negative
   set** and it's optional.

3. **Build the model files**. The minimum required parameters are:
   
   ```bash
   $> ROCker build -P positive.txt -o prep
   ```
   
   Where `positive.txt` is the set from step 1, and `prep` is the base name for the output files. You can also
   pass the negative set from step 2 using `-N` (or `-n`). For additional supported options, execute `ROCker build -h`.
   This is by far the most computationally-expensive step, so you might want to consider using multiple threads (`-t`)
   or even re-using files in case the run fails (`--reuse-files` and`--nocleanup`). Also, consider setting the
   simulated read length to match that of your metagenomes (`-l`).

4. **Compile the model**. The minimum required parameters are:
   
   ```bash
   $> ROCker compile -a prep.aln -b prep.blast -k model.rocker
   ```

   Where `prep.aln` is the alignment generated in step 3 (manual curation is strongly encouraged), `prep.blast` is the
   reference BLAST generated in step 3, and `model.rocker` is the model to compile.

5. **Register your model** (optional). If you would like to share your model with the community, please
   [Contact us](mailto:lrr@gatech.edu). We'll need the final ROCker model and the reference BLAST, and will add your
   model to our [curated list](http://enve-omics.ce.gatech.edu/rocker/models).



