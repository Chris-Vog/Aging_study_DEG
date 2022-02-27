

# Project: Extrinsic and intrinsic skin aging -Whole exon sequencing

Probanden wurden in drei Altersklassen (jung, mittel, alt) eingeteilt und biopsiert. Dabei wurden Proben von der Schulter (extrinsic) und dem Gesäß (intrinsic) entnommen. Die RNA wurde isoliert und in cDNA umgeschrieben. *Library preparation* und Sequenzierung mittels Illumina&reg; wurde extern durchgeführt. Die Indizes wurden in den bereitgestellten Proben bereits entfernt.

Die `URLs` wurden in der Textdatei `raw-data.txt` gespeichert und `wget` überführt.

  ```sh
  wget -P /path/to/Raw/ -i /path/to/txt-file
  ```

  ## Vorbereitung der Arbeitsumgebung mittels Anaconda

Anaconda ist eine Sammlung von Programmen, welches das Management und die Entwicklung von `Packages` vereinfacht und dabei hilft Abhängigkeiten von Programmen einzuhalten und zu verwalten. **Bioconda** ist dabei ein `Channel` für den **conda package manager**, der sich auf bioinformatische Software spezialisiert hat und die Installation dieser `Tools` stark vereinfacht. Zu Begin muss eine `Env` erstellt werden.

  ```bash
  conda create --name RNA-Seq #Der Name kann beliebig ausgewählt werden
  ```

Aktivierung der Arbeitsumgebung

  ```bash
  source activate RNA-Seq
  ```

Sobald die Arbeitsumgebung aktiviert ist, können die benötigten Programme innerhalb dieser Umgebung installiert werden. Der Vorteil bei der Verwendung von Anaconda ist, dass bei der Installation alle nötigen Abhängigkeiten ebenfalls installiert werden.

Diese `Env` kann darüber hinaus in einer `.yaml`-Datei gespeichert werden.

  ```bash
  conda env export --name env_name > ~/path/to/file.yaml #env_name mit eigentlichen Namen austauschen
  ```

Für die Reproduzierbarkeit der Ergebnisse, kann die so erstellte `.yaml`-Datei dazu verwendet werden die Umgebung auf weitere Maschinen zu klonen.

  ```bash
  conda env create -f /path/to/file.yml #Der Name der Umgebung wird automatisch übernommen
  ```

Die Installation erfolgt mit Hilfe des Befehls `conda install`:

  ```bash
  conda install -c bioconda fastqc #Mit -c wird der Channel angegeben, in dem sich das Programm befindet.
  ```

## Qualitätskontrolle

Die Qualität der Proben wurden mit Hilfe des Programms `fastqc` analysiert - ein Tool zur Analyse von High-Throughput Sequenzier-Daten.

  ```bash
  fastqc /path/to/fastq-files/*.fastq --outdir /path/to/folder
  ```

Die Ergebnisse werden als `.html`-Datei gespeichert. `*` wird hier als *wildcard* verwendet, sodass der Vorgang auf alle `.fastq`-Dateien angewendet wird.

Einige Daten konnten nicht prozessiert werden, da diese trunkiert waren. Der Fehler bei diesen Daten war darauf zurückzuführen, dass die Datenstruktur der `fastq`-Daten nicht gewahrt wurde.

## Data integrity

Auf dem Server befanden sich neben den `fastq`-Dateien auch `md5`-Checksums. Ich habe die Checksums der heruntergeladenen Dateien mit den bereitgestellten Checksums verglichen:

  ```bash
  md5sum /path/to/fastq
  ```

Die `checksums` waren identisch, was bedeutet, dass die Daten auf dem Server bereits fehlerhaft sind.

## Alignment (STAR)

### Download des Referenz-Genoms und der Annotation

Die Sequenzier-Daten werden bei diesem Schritt mit dem humanen Genom abgeglichen und die Position mit der höchsten Übereinstimmung identifiziert. Das humane Genom wird als `FASTA` zur Verfügung gestellt.

  ```bash
  wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/GRCh38.p13.genome.fa.gz
  ```

Nukleotid-Sequenz der GRCh38.p13 Genomassemblierung aller Regionen, inklusive Referenz-Chromosomen, Korrekturen und Haplotypen. Die dazugehörige Annotation ist in einer `.gtf`-Datei gespeichert. Wichtig hierbei ist, dass die Versionen übereinstimmen.

  ```bash
  wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.chr_patch_hapl_scaff.annotation.gtf.gz
  ```

### Alignment mit STAR

Spliced Transcripts Alignment to a Reference (STAR) ist ein schneller Mapper, der es ermöglicht auch Splice-Varianten und Fusion-Reads zu detektieren. Das Programm wird folgendermaßen installiert:

  ```bash
  conda install -c bioconda star
  ```

Bevor die einzelnen Sequenzen aligned werden können, muss ein Index vom Referenz-Genom erstellt werden:

  ```bash
  STAR \
  --runMode genomeGenerate \ #Argument zur Erstellung des Indexes
  --genomeDir genome/STAR/ \ #Ordner in den der Index gespeichert wird
  --genomeFastaFiles referenceGenome/GRCh38.p13.genome.fa \ #Pfad zum Referenzgenonm
  --sjdbGTFfile referenceGenome/gencode.v38.chr_patch_hapl_scaff.basic.annotation.gtf \ #Pfad zur Annotation
  --runThreadN 16
  
  Anschließend können die Sequenzen `aligned` werden:

  ```bash
  STAR \
  --genomeDir genome/STAR/ \ #Ordner mit indiziertem Genom
  --readFilesIn /path/to/fastq-Datei \
  --outFileNamePrefix /path/to/result-directory \
  --outSAMattributes All \ #Standard
  --outSAMunmapped Within \ #Wie mit nicht gemappten Reads verfahren
  --runThreadN 16 \ #Anzahl der zu verwendenen Threads
  --outSAMtype BAM SortedByCoordinate \ #In welchem Dateityp die Ergebnisse zu speichern sind und wie sie zu sortieren sind
  --readFilesCommand zcat #Welches Programm für verwendet wird um stdout zu generieren
  ```

Damit die Daten in Tools wie den `Integrative Genomic Viewer` visualisiert werden können, müssen die prozessierten Daten indiziert werden. Im Falle des `IGV` ist es notwendig, dass sich der Index und die `.bam`-Datei im selben Ordner befindet.

  ```bash
  samtools index /path/to/bam_file
  ```

Um alle `.bam`-Files automatisch zu indizieren, wurde folgender Befehl verwendet:

  ```bash
  #Der Pfad darf dabei nur zur höheren Ordner-Ebene führen
  find /path/to/bam_files -mindepth 1 -type d -exec bash -c "cd '{}' && samtools index Aligned.sortedByCoord.out.bam" \; 
  ```

### Qualitätskontrolle *Mapping*

Die Qualität und die Effizienz des *Mappings* wird ebenfalls mittels `samtools` analysiert. Hierzu werden die flagstats erstellt, die die Anzahl erfolgreicher Alignments aufzählen.

  ```bash
  samtools flagstat /path/to/bam_file > flagstat.txt
  ```

Die so erstellten `.txt`-Dateien werden mit Hilfe von **R** ausgewertet und visualisiert.

