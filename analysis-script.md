    knitr::opts_chunk$set(cache = TRUE, echo = FALSE)
    library(tidyverse)
    library(here)
    library(gt)
    library(rlang)
    library(minfi)

    select <- dplyr::select
    rename <- dplyr::rename
    filter <- dplyr::filter

    dir.create(here("figures"))

    set.seed(1234)


The **methcon5** package have been updated since this writeup. Please install old version to complete analysis


    #remotes::install_github("EmilHvitfeldt/methcon5@0.0.0.9000")
    library(methcon5)


Data Preparation
----------------

The following function, will take a folder containing .idat files and
calculates the noob adjusted beta values for the samples.

    idat_to_beta_values_noob <- function(base) {
      RGset <- read.metharray.exp(base = base)
      RGset@annotation = c(array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b4.hg19")

      MSet.noob <- preprocessNoob(RGset, offset = 15, dyeCorr = TRUE, verbose = TRUE)

      ratioSet.noob <- ratioConvert(MSet.noob, what =  "both", keepCN = TRUE)
      beta.noob <- getBeta(ratioSet.noob)
      beta.noob
    }

This function is applied to each of the tissue folders.

    colon_idat <- idat_to_beta_values_noob(base = "data-raw/colon-idat/")
    si_idat <- idat_to_beta_values_noob(base = "data-raw/si-idat/")
    endo_idat <- idat_to_beta_values_noob(base = "data-raw/endo-idat/")

Infinium MethylationEPIC v1.0 B4 Manifest File (CSV Format) is to be
downloaded from the following link
<a href="https://support.illumina.com/downloads/infinium-methylationepic-v1-0-product-files.html" class="uri">https://support.illumina.com/downloads/infinium-methylationepic-v1-0-product-files.html</a>
and `MethylationEPIC_v-1-0_B4.csv` should be placed in the folder
`data-raw`. `pwd()` will take a vector and calculate the average
manhattan distance. We then apply it row-wise to the methylation data.
Lastly we take the calulated average manhattan distances and append it
to the manifest file.

    pwd <- function(x) {
      mean(dist(x, method = "manhattan"))
    }

    colon_pwd_values <- colon_idat %>% apply(1, pwd)
    si_pwd_values <- si_idat %>% apply(1, pwd)
    endo_pwd_values <- endo_idat %>% apply(1, pwd)

    prepped_data <- 
      read_csv(here("data-raw/MethylationEPIC_v-1-0_B4.csv"), skip = 7, 
                         col_types = cols(Name = col_character(),
                                          CHR = col_character(),
                                          MAPINFO = col_double(),
                                          UCSC_RefGene_Name = col_character(),
                                          UCSC_RefGene_Group = col_character(),
                                          Relation_to_UCSC_CpG_Island = col_character())) %>%
      left_join(by = "Name",
                tibble(Name = names(colon_pwd_values),
                       colon_pwd = colon_pwd_values)) %>%
      left_join(by = "Name",
                tibble(Name = names(si_pwd_values),
                       si_pwd = si_pwd_values)) %>%
      left_join(by = "Name",
                tibble(Name = names(endo_pwd_values),
                       endo_pwd = endo_pwd_values)) %>%
      select(-IlmnID, -(AddressA_ID:Genome_Build), -(SourceSeq:Strand), -UCSC_RefGene_Accession, 
             -UCSC_CpG_Islands_Name, -Phantom4_Enhancers, -(DMR:X48), -Phantom5_Enhancers)

    write_csv(prepped_data, here("data", "prepped_data.csv"))

Next we seperate the rows to apply for a CpG site to appear in multiple
genes as annotated.

    gened_data <- prepped_data %>%
      arrange(CHR, MAPINFO) %>%
      mutate(temp = floor(row_number() / n() * 100)) %>%
      separate_rows(UCSC_RefGene_Name) %>% 
      distinct()

    write_csv(gened_data, here("data", "gened_data.csv"))

Here we apply the the first function from methcon5, we calculates the
unadjusted conservation for each gene for each of the tissue types.

    gened_data <- read_csv(here("data", "gened_data.csv"))
    pwd_by_gene <- gened_data %>%
      filter(!is.na(UCSC_RefGene_Name)) %>%
      ii_summarize(UCSC_RefGene_Name, colon_pwd:endo_pwd)

Next we apply the second function from methcon5 where we will apply two
different boot-strapping method to each of the tissuetypes.

    permed_data <- pwd_by_gene %>%
      filter(n < quantile(n, 0.99)) %>%
      bootter(gened_data, colon_pwd, UCSC_RefGene_Name, 1000,  method = 1) %>%
      bootter(gened_data, si_pwd,    UCSC_RefGene_Name, 1000, method = 1) %>%
      bootter(gened_data, endo_pwd,  UCSC_RefGene_Name, 1000, method = 1) %>%
      bootter(gened_data, colon_pwd, UCSC_RefGene_Name, 1000,  method = 3) %>%
      bootter(gened_data, si_pwd,    UCSC_RefGene_Name, 1000, method = 3) %>%
      bootter(gened_data, endo_pwd,  UCSC_RefGene_Name, 1000, method = 3)

    write_csv(permed_data, here("data", "permed_data.csv"))

    ## Parsed with column specification:
    ## cols(
    ##   Name = col_character(),
    ##   CHR = col_double(),
    ##   MAPINFO = col_double(),
    ##   UCSC_RefGene_Name = col_character(),
    ##   UCSC_RefGene_Group = col_character(),
    ##   Relation_to_UCSC_CpG_Island = col_character(),
    ##   colon_pwd = col_double(),
    ##   si_pwd = col_double(),
    ##   endo_pwd = col_double(),
    ##   temp = col_double()
    ## )

    ## Warning: 21500 parsing failures.
    ##    row col expected actual                                                                                                                  file
    ## 920898 CHR a double      X '/Users/emilhvitfeldthansen/Documents/USC/papers/Epigenetic-Conservation-Is-A-Beacon-Of-Function/data/gened_data.csv'
    ## 920899 CHR a double      X '/Users/emilhvitfeldthansen/Documents/USC/papers/Epigenetic-Conservation-Is-A-Beacon-Of-Function/data/gened_data.csv'
    ## 920900 CHR a double      X '/Users/emilhvitfeldthansen/Documents/USC/papers/Epigenetic-Conservation-Is-A-Beacon-Of-Function/data/gened_data.csv'
    ## 920901 CHR a double      X '/Users/emilhvitfeldthansen/Documents/USC/papers/Epigenetic-Conservation-Is-A-Beacon-Of-Function/data/gened_data.csv'
    ## 920902 CHR a double      X '/Users/emilhvitfeldthansen/Documents/USC/papers/Epigenetic-Conservation-Is-A-Beacon-Of-Function/data/gened_data.csv'
    ## ...... ... ........ ...... .....................................................................................................................
    ## See problems(...) for more details.

    ## Parsed with column specification:
    ## cols(
    ##   UCSC_RefGene_Name = col_character(),
    ##   colon_pwd = col_double(),
    ##   si_pwd = col_double(),
    ##   endo_pwd = col_double(),
    ##   n = col_double(),
    ##   colon_pwd_v1 = col_double(),
    ##   si_pwd_v1 = col_double(),
    ##   endo_pwd_v1 = col_double(),
    ##   colon_pwd_v3 = col_double(),
    ##   si_pwd_v3 = col_double(),
    ##   endo_pwd_v3 = col_double()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   Name = col_character(),
    ##   CHR = col_character(),
    ##   MAPINFO = col_double(),
    ##   UCSC_RefGene_Name = col_character(),
    ##   UCSC_RefGene_Group = col_character(),
    ##   Relation_to_UCSC_CpG_Island = col_character(),
    ##   colon_pwd = col_double(),
    ##   si_pwd = col_double(),
    ##   endo_pwd = col_double()
    ## )

Next we calculate the mean methylation for CpG sites in the promotor
region.

    promoted_data <- prepped_data %>%
      mutate(gene = str_extract(UCSC_RefGene_Name, "^[^;]+")) %>%
      filter(!is.na(gene)) %>%
      left_join(read_csv(here("data-raw/MethylationEPIC_v-1-0_B4.csv"), skip = 7) %>%  
                  dplyr::select(Name, Regulatory_Feature_Group), 
                by = "Name") %>%
      filter(Regulatory_Feature_Group == "Promoter_Associated") %>%
      group_by(gene) %>%
      summarise(colon_pro = mean(colon_pwd),
                si_pro = mean(si_pwd),
                endo_pro = mean(endo_pwd))

    ## Warning: Missing column names filled in: 'X48' [48]

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_character(),
    ##   Genome_Build = col_double(),
    ##   MAPINFO = col_double(),
    ##   `450k_Enhancer` = col_logical(),
    ##   DNase_Hypersensitivity_Evidence_Count = col_double(),
    ##   OpenChromatin_Evidence_Count = col_double(),
    ##   TFBS_Evidence_Count = col_double(),
    ##   Methyl27_Loci = col_logical(),
    ##   Methyl450_Loci = col_logical(),
    ##   Coordinate_36 = col_double(),
    ##   Random_Loci = col_logical()
    ## )

    ## See spec(...) for full column specifications.

    ## Warning: 866669 parsing failures.
    ## row col   expected     actual                                                                                                                                    file
    ##   1  -- 48 columns 47 columns '/Users/emilhvitfeldthansen/Documents/USC/papers/Epigenetic-Conservation-Is-A-Beacon-Of-Function/data-raw/MethylationEPIC_v-1-0_B4.csv'
    ##   2  -- 48 columns 47 columns '/Users/emilhvitfeldthansen/Documents/USC/papers/Epigenetic-Conservation-Is-A-Beacon-Of-Function/data-raw/MethylationEPIC_v-1-0_B4.csv'
    ##   3  -- 48 columns 47 columns '/Users/emilhvitfeldthansen/Documents/USC/papers/Epigenetic-Conservation-Is-A-Beacon-Of-Function/data-raw/MethylationEPIC_v-1-0_B4.csv'
    ##   4  -- 48 columns 47 columns '/Users/emilhvitfeldthansen/Documents/USC/papers/Epigenetic-Conservation-Is-A-Beacon-Of-Function/data-raw/MethylationEPIC_v-1-0_B4.csv'
    ##   5  -- 48 columns 47 columns '/Users/emilhvitfeldthansen/Documents/USC/papers/Epigenetic-Conservation-Is-A-Beacon-Of-Function/data-raw/MethylationEPIC_v-1-0_B4.csv'
    ## ... ... .......... .......... .......................................................................................................................................
    ## See problems(...) for more details.

Expression data can be downloaded from the following website
<a href="https://www.ebi.ac.uk/gxa/experiments/E-MTAB-2836/Results" class="uri">https://www.ebi.ac.uk/gxa/experiments/E-MTAB-2836/Results</a>,
and the resulting file `E-MTAB-2836-query-results.tpms.tsv` should be
placed in the folder named `data-raw`.

    expression_atlas <- read_tsv(here("data-raw", "E-MTAB-2836-query-results.tpms.tsv"), skip = 4, 
                                 col_types = cols(
      .default = col_double(),
      `Gene ID` = col_character(),
      `Gene Name` = col_character()
    ))

Gene enrichment set analysis
----------------------------

In the following code we apply gene enrichment set analysis on the 5%
most conserved genes for each tissue type one at a time.

    tidy_gsea <- function(genes) {
      AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, genes, 'ENTREZID', 'SYMBOL') %>%
      ReactomePA::enrichPathway()
    }

    enrich_colon <- permed_data %>%
      filter(colon_pwd_v3 > quantile(colon_pwd_v3, 0.95, na.rm = TRUE)) %>%
      pull(UCSC_RefGene_Name) %>%
      tidy_gsea()

    ## 

    ## 'select()' returned 1:1 mapping between keys and columns

    ## 

    ## Registered S3 method overwritten by 'enrichplot':
    ##   method               from
    ##   fortify.enrichResult DOSE

    ## Loading required package: org.Hs.eg.db

    ## Loading required package: AnnotationDbi

    ## Warning: package 'AnnotationDbi' was built under R version 3.6.1

    ## 
    ## Attaching package: 'AnnotationDbi'

    ## The following object is masked _by_ '.GlobalEnv':
    ## 
    ##     select

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     select

    enrich_si <- permed_data %>%
      filter(si_pwd_v3 > quantile(si_pwd_v3, 0.95, na.rm = TRUE)) %>%
      pull(UCSC_RefGene_Name) %>%
      tidy_gsea()

    ## 'select()' returned 1:many mapping between keys and columns

    enrich_endo <- permed_data %>%
      filter(endo_pwd_v3 > quantile(endo_pwd_v3, 0.95, na.rm = TRUE)) %>%
      pull(UCSC_RefGene_Name) %>%
      tidy_gsea()

    ## 'select()' returned 1:1 mapping between keys and columns

Table 1
-------

The calculations happening in the following bits of code follow the same
steps. Split the data accoding to a specified region. Then calculate the
mean pwd for each tissue type.

    prepped_data %>%
      group_by(`Gene Associated` = !is.na(UCSC_RefGene_Name)) %>%
      summarise_at(vars(colon_pwd:endo_pwd), mean, na.rm = TRUE) %>%
      rename(Colon = colon_pwd,
             `Small Intestine` = si_pwd,
             `Endometrium` = endo_pwd) %>%
      mutate(`Gene Associated` = factor(`Gene Associated`, c(TRUE, FALSE), c("Yes", "No"))) %>%
      arrange(rev(rownames(.))) %>%
      gt() %>%
      fmt_number(columns = vars(Colon, `Small Intestine`, `Endometrium`), decimals = 3) %>%
      tab_spanner(
        label = "Manhattan distance",
        columns = vars(Colon, `Small Intestine`, `Endometrium`)
      )

<!--html_preserve-->
<style>html {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
}

#pzveloiefp .gt_table {
  display: table;
  border-collapse: collapse;
  margin-left: auto;
  margin-right: auto;
  color: #000000;
  font-size: 16px;
  background-color: #FFFFFF;
  /* table.background.color */
  width: auto;
  /* table.width */
  border-top-style: solid;
  /* table.border.top.style */
  border-top-width: 2px;
  /* table.border.top.width */
  border-top-color: #A8A8A8;
  /* table.border.top.color */
}

#pzveloiefp .gt_heading {
  background-color: #FFFFFF;
  /* heading.background.color */
  border-bottom-color: #FFFFFF;
}

#pzveloiefp .gt_title {
  color: #000000;
  font-size: 125%;
  /* heading.title.font.size */
  padding-top: 4px;
  /* heading.top.padding */
  padding-bottom: 1px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#pzveloiefp .gt_subtitle {
  color: #000000;
  font-size: 85%;
  /* heading.subtitle.font.size */
  padding-top: 1px;
  padding-bottom: 4px;
  /* heading.bottom.padding */
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#pzveloiefp .gt_bottom_border {
  border-bottom-style: solid;
  /* heading.border.bottom.style */
  border-bottom-width: 2px;
  /* heading.border.bottom.width */
  border-bottom-color: #A8A8A8;
  /* heading.border.bottom.color */
}

#pzveloiefp .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  padding-top: 4px;
  padding-bottom: 4px;
}

#pzveloiefp .gt_col_heading {
  color: #000000;
  background-color: #FFFFFF;
  /* column_labels.background.color */
  font-size: 16px;
  /* column_labels.font.size */
  font-weight: initial;
  /* column_labels.font.weight */
  vertical-align: middle;
  padding: 10px;
  margin: 10px;
}

#pzveloiefp .gt_sep_right {
  border-right: 5px solid #FFFFFF;
}

#pzveloiefp .gt_group_heading {
  padding: 8px;
  color: #000000;
  background-color: #FFFFFF;
  /* row_group.background.color */
  font-size: 16px;
  /* row_group.font.size */
  font-weight: initial;
  /* row_group.font.weight */
  border-top-style: solid;
  /* row_group.border.top.style */
  border-top-width: 2px;
  /* row_group.border.top.width */
  border-top-color: #A8A8A8;
  /* row_group.border.top.color */
  border-bottom-style: solid;
  /* row_group.border.bottom.style */
  border-bottom-width: 2px;
  /* row_group.border.bottom.width */
  border-bottom-color: #A8A8A8;
  /* row_group.border.bottom.color */
  vertical-align: middle;
}

#pzveloiefp .gt_empty_group_heading {
  padding: 0.5px;
  color: #000000;
  background-color: #FFFFFF;
  /* row_group.background.color */
  font-size: 16px;
  /* row_group.font.size */
  font-weight: initial;
  /* row_group.font.weight */
  border-top-style: solid;
  /* row_group.border.top.style */
  border-top-width: 2px;
  /* row_group.border.top.width */
  border-top-color: #A8A8A8;
  /* row_group.border.top.color */
  border-bottom-style: solid;
  /* row_group.border.bottom.style */
  border-bottom-width: 2px;
  /* row_group.border.bottom.width */
  border-bottom-color: #A8A8A8;
  /* row_group.border.bottom.color */
  vertical-align: middle;
}

#pzveloiefp .gt_striped {
  background-color: #f2f2f2;
}

#pzveloiefp .gt_from_md > :first-child {
  margin-top: 0;
}

#pzveloiefp .gt_from_md > :last-child {
  margin-bottom: 0;
}

#pzveloiefp .gt_row {
  padding: 8px;
  /* row.padding */
  margin: 10px;
  vertical-align: middle;
}

#pzveloiefp .gt_stub {
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #A8A8A8;
  padding-left: 12px;
}

#pzveloiefp .gt_summary_row {
  color: #000000;
  background-color: #FFFFFF;
  /* summary_row.background.color */
  padding: 8px;
  /* summary_row.padding */
  text-transform: inherit;
  /* summary_row.text_transform */
}

#pzveloiefp .gt_grand_summary_row {
  color: #000000;
  background-color: #FFFFFF;
  /* grand_summary_row.background.color */
  padding: 8px;
  /* grand_summary_row.padding */
  text-transform: inherit;
  /* grand_summary_row.text_transform */
}

#pzveloiefp .gt_first_summary_row {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
}

#pzveloiefp .gt_first_grand_summary_row {
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #A8A8A8;
}

#pzveloiefp .gt_table_body {
  border-top-style: solid;
  /* table_body.border.top.style */
  border-top-width: 2px;
  /* table_body.border.top.width */
  border-top-color: #A8A8A8;
  /* table_body.border.top.color */
  border-bottom-style: solid;
  /* table_body.border.bottom.style */
  border-bottom-width: 2px;
  /* table_body.border.bottom.width */
  border-bottom-color: #A8A8A8;
  /* table_body.border.bottom.color */
}

#pzveloiefp .gt_footnote {
  font-size: 90%;
  /* footnote.font.size */
  padding: 4px;
  /* footnote.padding */
}

#pzveloiefp .gt_sourcenote {
  font-size: 90%;
  /* sourcenote.font.size */
  padding: 4px;
  /* sourcenote.padding */
}

#pzveloiefp .gt_center {
  text-align: center;
}

#pzveloiefp .gt_left {
  text-align: left;
}

#pzveloiefp .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#pzveloiefp .gt_font_normal {
  font-weight: normal;
}

#pzveloiefp .gt_font_bold {
  font-weight: bold;
}

#pzveloiefp .gt_font_italic {
  font-style: italic;
}

#pzveloiefp .gt_super {
  font-size: 65%;
}

#pzveloiefp .gt_footnote_glyph {
  font-style: italic;
  font-size: 65%;
}
</style>

<!--gt table start-->
<table class="gt_table">
<tr>
<th class="gt_col_heading gt_center" rowspan="2" colspan="1">
Gene Associated
</th>
<th class="gt_col_heading gt_column_spanner gt_center" rowspan="1" colspan="3">
Manhattan distance
</th>
</tr>
<tr>
<th class="gt_col_heading gt_right" rowspan="1" colspan="1">
Colon
</th>
<th class="gt_col_heading gt_right" rowspan="1" colspan="1">
Small Intestine
</th>
<th class="gt_col_heading gt_NA" rowspan="1" colspan="1">
Endometrium
</th>
</tr>
<tbody class="gt_table_body">
<tr>
<td class="gt_row gt_center">
Yes
</td>
<td class="gt_row gt_right">
0.088
</td>
<td class="gt_row gt_right">
0.101
</td>
<td class="gt_row gt_right">
0.081
</td>
</tr>
<tr>
<td class="gt_row gt_center gt_striped">
No
</td>
<td class="gt_row gt_right gt_striped">
0.118
</td>
<td class="gt_row gt_right gt_striped">
0.134
</td>
<td class="gt_row gt_right gt_striped">
0.112
</td>
</tr>
</tbody>
</table>
<!--gt table end-->

<!--/html_preserve-->

    prepped_data %>%
      group_by(`Island relation` = Relation_to_UCSC_CpG_Island) %>%
      summarise_at(vars(colon_pwd:endo_pwd), mean, na.rm = TRUE) %>%
      rename(Colon = colon_pwd,
             `Small Intestine` = si_pwd,
             `Endometrium` = endo_pwd) %>%
      dplyr::rename(ir = `Island relation`) %>%
      mutate(ir = case_when(is.na(ir) ~ "Sea",
                            str_detect(ir, "^S") ~ str_replace(ir, "^S_", "South "),
                            str_detect(ir, "^N") ~ str_replace(ir, "^N_", "North "),
                            TRUE ~ ir)) %>%
      dplyr::rename(`Island relation` = ir) %>%
      arrange(rev(rownames(.))) %>%
      gt() %>%
      fmt_number(columns = vars(Colon, `Small Intestine`, `Endometrium`), decimals = 3) %>%
      tab_spanner(
        label = "Manhattan distance",
        columns = vars(Colon, `Small Intestine`, `Endometrium`)
      )

<!--html_preserve-->
<style>html {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
}

#dbgvzfontn .gt_table {
  display: table;
  border-collapse: collapse;
  margin-left: auto;
  margin-right: auto;
  color: #000000;
  font-size: 16px;
  background-color: #FFFFFF;
  /* table.background.color */
  width: auto;
  /* table.width */
  border-top-style: solid;
  /* table.border.top.style */
  border-top-width: 2px;
  /* table.border.top.width */
  border-top-color: #A8A8A8;
  /* table.border.top.color */
}

#dbgvzfontn .gt_heading {
  background-color: #FFFFFF;
  /* heading.background.color */
  border-bottom-color: #FFFFFF;
}

#dbgvzfontn .gt_title {
  color: #000000;
  font-size: 125%;
  /* heading.title.font.size */
  padding-top: 4px;
  /* heading.top.padding */
  padding-bottom: 1px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#dbgvzfontn .gt_subtitle {
  color: #000000;
  font-size: 85%;
  /* heading.subtitle.font.size */
  padding-top: 1px;
  padding-bottom: 4px;
  /* heading.bottom.padding */
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#dbgvzfontn .gt_bottom_border {
  border-bottom-style: solid;
  /* heading.border.bottom.style */
  border-bottom-width: 2px;
  /* heading.border.bottom.width */
  border-bottom-color: #A8A8A8;
  /* heading.border.bottom.color */
}

#dbgvzfontn .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  padding-top: 4px;
  padding-bottom: 4px;
}

#dbgvzfontn .gt_col_heading {
  color: #000000;
  background-color: #FFFFFF;
  /* column_labels.background.color */
  font-size: 16px;
  /* column_labels.font.size */
  font-weight: initial;
  /* column_labels.font.weight */
  vertical-align: middle;
  padding: 10px;
  margin: 10px;
}

#dbgvzfontn .gt_sep_right {
  border-right: 5px solid #FFFFFF;
}

#dbgvzfontn .gt_group_heading {
  padding: 8px;
  color: #000000;
  background-color: #FFFFFF;
  /* row_group.background.color */
  font-size: 16px;
  /* row_group.font.size */
  font-weight: initial;
  /* row_group.font.weight */
  border-top-style: solid;
  /* row_group.border.top.style */
  border-top-width: 2px;
  /* row_group.border.top.width */
  border-top-color: #A8A8A8;
  /* row_group.border.top.color */
  border-bottom-style: solid;
  /* row_group.border.bottom.style */
  border-bottom-width: 2px;
  /* row_group.border.bottom.width */
  border-bottom-color: #A8A8A8;
  /* row_group.border.bottom.color */
  vertical-align: middle;
}

#dbgvzfontn .gt_empty_group_heading {
  padding: 0.5px;
  color: #000000;
  background-color: #FFFFFF;
  /* row_group.background.color */
  font-size: 16px;
  /* row_group.font.size */
  font-weight: initial;
  /* row_group.font.weight */
  border-top-style: solid;
  /* row_group.border.top.style */
  border-top-width: 2px;
  /* row_group.border.top.width */
  border-top-color: #A8A8A8;
  /* row_group.border.top.color */
  border-bottom-style: solid;
  /* row_group.border.bottom.style */
  border-bottom-width: 2px;
  /* row_group.border.bottom.width */
  border-bottom-color: #A8A8A8;
  /* row_group.border.bottom.color */
  vertical-align: middle;
}

#dbgvzfontn .gt_striped {
  background-color: #f2f2f2;
}

#dbgvzfontn .gt_from_md > :first-child {
  margin-top: 0;
}

#dbgvzfontn .gt_from_md > :last-child {
  margin-bottom: 0;
}

#dbgvzfontn .gt_row {
  padding: 8px;
  /* row.padding */
  margin: 10px;
  vertical-align: middle;
}

#dbgvzfontn .gt_stub {
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #A8A8A8;
  padding-left: 12px;
}

#dbgvzfontn .gt_summary_row {
  color: #000000;
  background-color: #FFFFFF;
  /* summary_row.background.color */
  padding: 8px;
  /* summary_row.padding */
  text-transform: inherit;
  /* summary_row.text_transform */
}

#dbgvzfontn .gt_grand_summary_row {
  color: #000000;
  background-color: #FFFFFF;
  /* grand_summary_row.background.color */
  padding: 8px;
  /* grand_summary_row.padding */
  text-transform: inherit;
  /* grand_summary_row.text_transform */
}

#dbgvzfontn .gt_first_summary_row {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
}

#dbgvzfontn .gt_first_grand_summary_row {
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #A8A8A8;
}

#dbgvzfontn .gt_table_body {
  border-top-style: solid;
  /* table_body.border.top.style */
  border-top-width: 2px;
  /* table_body.border.top.width */
  border-top-color: #A8A8A8;
  /* table_body.border.top.color */
  border-bottom-style: solid;
  /* table_body.border.bottom.style */
  border-bottom-width: 2px;
  /* table_body.border.bottom.width */
  border-bottom-color: #A8A8A8;
  /* table_body.border.bottom.color */
}

#dbgvzfontn .gt_footnote {
  font-size: 90%;
  /* footnote.font.size */
  padding: 4px;
  /* footnote.padding */
}

#dbgvzfontn .gt_sourcenote {
  font-size: 90%;
  /* sourcenote.font.size */
  padding: 4px;
  /* sourcenote.padding */
}

#dbgvzfontn .gt_center {
  text-align: center;
}

#dbgvzfontn .gt_left {
  text-align: left;
}

#dbgvzfontn .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#dbgvzfontn .gt_font_normal {
  font-weight: normal;
}

#dbgvzfontn .gt_font_bold {
  font-weight: bold;
}

#dbgvzfontn .gt_font_italic {
  font-style: italic;
}

#dbgvzfontn .gt_super {
  font-size: 65%;
}

#dbgvzfontn .gt_footnote_glyph {
  font-style: italic;
  font-size: 65%;
}
</style>

<!--gt table start-->
<table class="gt_table">
<tr>
<th class="gt_col_heading gt_center" rowspan="2" colspan="1">
Island relation
</th>
<th class="gt_col_heading gt_column_spanner gt_center" rowspan="1" colspan="3">
Manhattan distance
</th>
</tr>
<tr>
<th class="gt_col_heading gt_right" rowspan="1" colspan="1">
Colon
</th>
<th class="gt_col_heading gt_right" rowspan="1" colspan="1">
Small Intestine
</th>
<th class="gt_col_heading gt_NA" rowspan="1" colspan="1">
Endometrium
</th>
</tr>
<tbody class="gt_table_body">
<tr>
<td class="gt_row gt_left">
Sea
</td>
<td class="gt_row gt_right">
0.107
</td>
<td class="gt_row gt_right">
0.125
</td>
<td class="gt_row gt_right">
0.102
</td>
</tr>
<tr>
<td class="gt_row gt_left gt_striped">
South Shore
</td>
<td class="gt_row gt_right gt_striped">
0.102
</td>
<td class="gt_row gt_right gt_striped">
0.112
</td>
<td class="gt_row gt_right gt_striped">
0.089
</td>
</tr>
<tr>
<td class="gt_row gt_left">
South Shelf
</td>
<td class="gt_row gt_right">
0.107
</td>
<td class="gt_row gt_right">
0.121
</td>
<td class="gt_row gt_right">
0.097
</td>
</tr>
<tr>
<td class="gt_row gt_left gt_striped">
North Shore
</td>
<td class="gt_row gt_right gt_striped">
0.105
</td>
<td class="gt_row gt_right gt_striped">
0.115
</td>
<td class="gt_row gt_right gt_striped">
0.092
</td>
</tr>
<tr>
<td class="gt_row gt_left">
North Shelf
</td>
<td class="gt_row gt_right">
0.107
</td>
<td class="gt_row gt_right">
0.121
</td>
<td class="gt_row gt_right">
0.097
</td>
</tr>
<tr>
<td class="gt_row gt_left gt_striped">
Island
</td>
<td class="gt_row gt_right gt_striped">
0.055
</td>
<td class="gt_row gt_right gt_striped">
0.060
</td>
<td class="gt_row gt_right gt_striped">
0.051
</td>
</tr>
</tbody>
</table>
<!--gt table end-->

<!--/html_preserve-->

    prepped_data %>%
      group_by(`5'UTR` = str_detect(UCSC_RefGene_Group, "5'UTR") & !is.na(UCSC_RefGene_Group)) %>%
      summarise_at(vars(colon_pwd:endo_pwd), mean, na.rm = TRUE) %>%
      rename(Colon = colon_pwd,
             `Small Intestine` = si_pwd,
             `Endometrium` = endo_pwd) %>%
      mutate(`5'UTR` = factor(`5'UTR`, c(TRUE, FALSE), c("Yes", "No"))) %>%
      arrange(rev(rownames(.))) %>%
      gt() %>%
      fmt_number(columns = vars(Colon, `Small Intestine`, `Endometrium`), decimals = 3) %>%
      tab_spanner(
        label = "Manhattan distance",
        columns = vars(Colon, `Small Intestine`, `Endometrium`)
      )

<!--html_preserve-->
<style>html {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
}

#xdduhtxcdz .gt_table {
  display: table;
  border-collapse: collapse;
  margin-left: auto;
  margin-right: auto;
  color: #000000;
  font-size: 16px;
  background-color: #FFFFFF;
  /* table.background.color */
  width: auto;
  /* table.width */
  border-top-style: solid;
  /* table.border.top.style */
  border-top-width: 2px;
  /* table.border.top.width */
  border-top-color: #A8A8A8;
  /* table.border.top.color */
}

#xdduhtxcdz .gt_heading {
  background-color: #FFFFFF;
  /* heading.background.color */
  border-bottom-color: #FFFFFF;
}

#xdduhtxcdz .gt_title {
  color: #000000;
  font-size: 125%;
  /* heading.title.font.size */
  padding-top: 4px;
  /* heading.top.padding */
  padding-bottom: 1px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#xdduhtxcdz .gt_subtitle {
  color: #000000;
  font-size: 85%;
  /* heading.subtitle.font.size */
  padding-top: 1px;
  padding-bottom: 4px;
  /* heading.bottom.padding */
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#xdduhtxcdz .gt_bottom_border {
  border-bottom-style: solid;
  /* heading.border.bottom.style */
  border-bottom-width: 2px;
  /* heading.border.bottom.width */
  border-bottom-color: #A8A8A8;
  /* heading.border.bottom.color */
}

#xdduhtxcdz .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  padding-top: 4px;
  padding-bottom: 4px;
}

#xdduhtxcdz .gt_col_heading {
  color: #000000;
  background-color: #FFFFFF;
  /* column_labels.background.color */
  font-size: 16px;
  /* column_labels.font.size */
  font-weight: initial;
  /* column_labels.font.weight */
  vertical-align: middle;
  padding: 10px;
  margin: 10px;
}

#xdduhtxcdz .gt_sep_right {
  border-right: 5px solid #FFFFFF;
}

#xdduhtxcdz .gt_group_heading {
  padding: 8px;
  color: #000000;
  background-color: #FFFFFF;
  /* row_group.background.color */
  font-size: 16px;
  /* row_group.font.size */
  font-weight: initial;
  /* row_group.font.weight */
  border-top-style: solid;
  /* row_group.border.top.style */
  border-top-width: 2px;
  /* row_group.border.top.width */
  border-top-color: #A8A8A8;
  /* row_group.border.top.color */
  border-bottom-style: solid;
  /* row_group.border.bottom.style */
  border-bottom-width: 2px;
  /* row_group.border.bottom.width */
  border-bottom-color: #A8A8A8;
  /* row_group.border.bottom.color */
  vertical-align: middle;
}

#xdduhtxcdz .gt_empty_group_heading {
  padding: 0.5px;
  color: #000000;
  background-color: #FFFFFF;
  /* row_group.background.color */
  font-size: 16px;
  /* row_group.font.size */
  font-weight: initial;
  /* row_group.font.weight */
  border-top-style: solid;
  /* row_group.border.top.style */
  border-top-width: 2px;
  /* row_group.border.top.width */
  border-top-color: #A8A8A8;
  /* row_group.border.top.color */
  border-bottom-style: solid;
  /* row_group.border.bottom.style */
  border-bottom-width: 2px;
  /* row_group.border.bottom.width */
  border-bottom-color: #A8A8A8;
  /* row_group.border.bottom.color */
  vertical-align: middle;
}

#xdduhtxcdz .gt_striped {
  background-color: #f2f2f2;
}

#xdduhtxcdz .gt_from_md > :first-child {
  margin-top: 0;
}

#xdduhtxcdz .gt_from_md > :last-child {
  margin-bottom: 0;
}

#xdduhtxcdz .gt_row {
  padding: 8px;
  /* row.padding */
  margin: 10px;
  vertical-align: middle;
}

#xdduhtxcdz .gt_stub {
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #A8A8A8;
  padding-left: 12px;
}

#xdduhtxcdz .gt_summary_row {
  color: #000000;
  background-color: #FFFFFF;
  /* summary_row.background.color */
  padding: 8px;
  /* summary_row.padding */
  text-transform: inherit;
  /* summary_row.text_transform */
}

#xdduhtxcdz .gt_grand_summary_row {
  color: #000000;
  background-color: #FFFFFF;
  /* grand_summary_row.background.color */
  padding: 8px;
  /* grand_summary_row.padding */
  text-transform: inherit;
  /* grand_summary_row.text_transform */
}

#xdduhtxcdz .gt_first_summary_row {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
}

#xdduhtxcdz .gt_first_grand_summary_row {
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #A8A8A8;
}

#xdduhtxcdz .gt_table_body {
  border-top-style: solid;
  /* table_body.border.top.style */
  border-top-width: 2px;
  /* table_body.border.top.width */
  border-top-color: #A8A8A8;
  /* table_body.border.top.color */
  border-bottom-style: solid;
  /* table_body.border.bottom.style */
  border-bottom-width: 2px;
  /* table_body.border.bottom.width */
  border-bottom-color: #A8A8A8;
  /* table_body.border.bottom.color */
}

#xdduhtxcdz .gt_footnote {
  font-size: 90%;
  /* footnote.font.size */
  padding: 4px;
  /* footnote.padding */
}

#xdduhtxcdz .gt_sourcenote {
  font-size: 90%;
  /* sourcenote.font.size */
  padding: 4px;
  /* sourcenote.padding */
}

#xdduhtxcdz .gt_center {
  text-align: center;
}

#xdduhtxcdz .gt_left {
  text-align: left;
}

#xdduhtxcdz .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#xdduhtxcdz .gt_font_normal {
  font-weight: normal;
}

#xdduhtxcdz .gt_font_bold {
  font-weight: bold;
}

#xdduhtxcdz .gt_font_italic {
  font-style: italic;
}

#xdduhtxcdz .gt_super {
  font-size: 65%;
}

#xdduhtxcdz .gt_footnote_glyph {
  font-style: italic;
  font-size: 65%;
}
</style>

<!--gt table start-->
<table class="gt_table">
<tr>
<th class="gt_col_heading gt_center" rowspan="2" colspan="1">
5’UTR
</th>
<th class="gt_col_heading gt_column_spanner gt_center" rowspan="1" colspan="3">
Manhattan distance
</th>
</tr>
<tr>
<th class="gt_col_heading gt_right" rowspan="1" colspan="1">
Colon
</th>
<th class="gt_col_heading gt_right" rowspan="1" colspan="1">
Small Intestine
</th>
<th class="gt_col_heading gt_NA" rowspan="1" colspan="1">
Endometrium
</th>
</tr>
<tbody class="gt_table_body">
<tr>
<td class="gt_row gt_center">
Yes
</td>
<td class="gt_row gt_right">
0.079
</td>
<td class="gt_row gt_right">
0.092
</td>
<td class="gt_row gt_right">
0.072
</td>
</tr>
<tr>
<td class="gt_row gt_center gt_striped">
No
</td>
<td class="gt_row gt_right gt_striped">
0.099
</td>
<td class="gt_row gt_right gt_striped">
0.113
</td>
<td class="gt_row gt_right gt_striped">
0.093
</td>
</tr>
</tbody>
</table>
<!--gt table end-->

<!--/html_preserve-->

    prepped_data %>%
      group_by(`TSS1500` = str_detect(UCSC_RefGene_Group, "TSS1500") & !is.na(UCSC_RefGene_Group)) %>%
      summarise_at(vars(colon_pwd:endo_pwd), mean, na.rm = TRUE) %>%
      rename(Colon = colon_pwd,
             `Small Intestine` = si_pwd,
             `Endometrium` = endo_pwd) %>%
      mutate(`TSS1500` = factor(`TSS1500`, c(TRUE, FALSE), c("Yes", "No"))) %>%
      arrange(rev(rownames(.))) %>%
      gt() %>%
      fmt_number(columns = vars(Colon, `Small Intestine`, `Endometrium`), decimals = 3) %>%
      tab_spanner(
        label = "Manhattan distance",
        columns = vars(Colon, `Small Intestine`, `Endometrium`)
      )

<!--html_preserve-->
<style>html {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
}

#ebohtplcwi .gt_table {
  display: table;
  border-collapse: collapse;
  margin-left: auto;
  margin-right: auto;
  color: #000000;
  font-size: 16px;
  background-color: #FFFFFF;
  /* table.background.color */
  width: auto;
  /* table.width */
  border-top-style: solid;
  /* table.border.top.style */
  border-top-width: 2px;
  /* table.border.top.width */
  border-top-color: #A8A8A8;
  /* table.border.top.color */
}

#ebohtplcwi .gt_heading {
  background-color: #FFFFFF;
  /* heading.background.color */
  border-bottom-color: #FFFFFF;
}

#ebohtplcwi .gt_title {
  color: #000000;
  font-size: 125%;
  /* heading.title.font.size */
  padding-top: 4px;
  /* heading.top.padding */
  padding-bottom: 1px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#ebohtplcwi .gt_subtitle {
  color: #000000;
  font-size: 85%;
  /* heading.subtitle.font.size */
  padding-top: 1px;
  padding-bottom: 4px;
  /* heading.bottom.padding */
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#ebohtplcwi .gt_bottom_border {
  border-bottom-style: solid;
  /* heading.border.bottom.style */
  border-bottom-width: 2px;
  /* heading.border.bottom.width */
  border-bottom-color: #A8A8A8;
  /* heading.border.bottom.color */
}

#ebohtplcwi .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  padding-top: 4px;
  padding-bottom: 4px;
}

#ebohtplcwi .gt_col_heading {
  color: #000000;
  background-color: #FFFFFF;
  /* column_labels.background.color */
  font-size: 16px;
  /* column_labels.font.size */
  font-weight: initial;
  /* column_labels.font.weight */
  vertical-align: middle;
  padding: 10px;
  margin: 10px;
}

#ebohtplcwi .gt_sep_right {
  border-right: 5px solid #FFFFFF;
}

#ebohtplcwi .gt_group_heading {
  padding: 8px;
  color: #000000;
  background-color: #FFFFFF;
  /* row_group.background.color */
  font-size: 16px;
  /* row_group.font.size */
  font-weight: initial;
  /* row_group.font.weight */
  border-top-style: solid;
  /* row_group.border.top.style */
  border-top-width: 2px;
  /* row_group.border.top.width */
  border-top-color: #A8A8A8;
  /* row_group.border.top.color */
  border-bottom-style: solid;
  /* row_group.border.bottom.style */
  border-bottom-width: 2px;
  /* row_group.border.bottom.width */
  border-bottom-color: #A8A8A8;
  /* row_group.border.bottom.color */
  vertical-align: middle;
}

#ebohtplcwi .gt_empty_group_heading {
  padding: 0.5px;
  color: #000000;
  background-color: #FFFFFF;
  /* row_group.background.color */
  font-size: 16px;
  /* row_group.font.size */
  font-weight: initial;
  /* row_group.font.weight */
  border-top-style: solid;
  /* row_group.border.top.style */
  border-top-width: 2px;
  /* row_group.border.top.width */
  border-top-color: #A8A8A8;
  /* row_group.border.top.color */
  border-bottom-style: solid;
  /* row_group.border.bottom.style */
  border-bottom-width: 2px;
  /* row_group.border.bottom.width */
  border-bottom-color: #A8A8A8;
  /* row_group.border.bottom.color */
  vertical-align: middle;
}

#ebohtplcwi .gt_striped {
  background-color: #f2f2f2;
}

#ebohtplcwi .gt_from_md > :first-child {
  margin-top: 0;
}

#ebohtplcwi .gt_from_md > :last-child {
  margin-bottom: 0;
}

#ebohtplcwi .gt_row {
  padding: 8px;
  /* row.padding */
  margin: 10px;
  vertical-align: middle;
}

#ebohtplcwi .gt_stub {
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #A8A8A8;
  padding-left: 12px;
}

#ebohtplcwi .gt_summary_row {
  color: #000000;
  background-color: #FFFFFF;
  /* summary_row.background.color */
  padding: 8px;
  /* summary_row.padding */
  text-transform: inherit;
  /* summary_row.text_transform */
}

#ebohtplcwi .gt_grand_summary_row {
  color: #000000;
  background-color: #FFFFFF;
  /* grand_summary_row.background.color */
  padding: 8px;
  /* grand_summary_row.padding */
  text-transform: inherit;
  /* grand_summary_row.text_transform */
}

#ebohtplcwi .gt_first_summary_row {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
}

#ebohtplcwi .gt_first_grand_summary_row {
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #A8A8A8;
}

#ebohtplcwi .gt_table_body {
  border-top-style: solid;
  /* table_body.border.top.style */
  border-top-width: 2px;
  /* table_body.border.top.width */
  border-top-color: #A8A8A8;
  /* table_body.border.top.color */
  border-bottom-style: solid;
  /* table_body.border.bottom.style */
  border-bottom-width: 2px;
  /* table_body.border.bottom.width */
  border-bottom-color: #A8A8A8;
  /* table_body.border.bottom.color */
}

#ebohtplcwi .gt_footnote {
  font-size: 90%;
  /* footnote.font.size */
  padding: 4px;
  /* footnote.padding */
}

#ebohtplcwi .gt_sourcenote {
  font-size: 90%;
  /* sourcenote.font.size */
  padding: 4px;
  /* sourcenote.padding */
}

#ebohtplcwi .gt_center {
  text-align: center;
}

#ebohtplcwi .gt_left {
  text-align: left;
}

#ebohtplcwi .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#ebohtplcwi .gt_font_normal {
  font-weight: normal;
}

#ebohtplcwi .gt_font_bold {
  font-weight: bold;
}

#ebohtplcwi .gt_font_italic {
  font-style: italic;
}

#ebohtplcwi .gt_super {
  font-size: 65%;
}

#ebohtplcwi .gt_footnote_glyph {
  font-style: italic;
  font-size: 65%;
}
</style>

<!--gt table start-->
<table class="gt_table">
<tr>
<th class="gt_col_heading gt_center" rowspan="2" colspan="1">
TSS1500
</th>
<th class="gt_col_heading gt_column_spanner gt_center" rowspan="1" colspan="3">
Manhattan distance
</th>
</tr>
<tr>
<th class="gt_col_heading gt_right" rowspan="1" colspan="1">
Colon
</th>
<th class="gt_col_heading gt_right" rowspan="1" colspan="1">
Small Intestine
</th>
<th class="gt_col_heading gt_NA" rowspan="1" colspan="1">
Endometrium
</th>
</tr>
<tbody class="gt_table_body">
<tr>
<td class="gt_row gt_center">
Yes
</td>
<td class="gt_row gt_right">
0.089
</td>
<td class="gt_row gt_right">
0.097
</td>
<td class="gt_row gt_right">
0.080
</td>
</tr>
<tr>
<td class="gt_row gt_center gt_striped">
No
</td>
<td class="gt_row gt_right gt_striped">
0.098
</td>
<td class="gt_row gt_right gt_striped">
0.113
</td>
<td class="gt_row gt_right gt_striped">
0.092
</td>
</tr>
</tbody>
</table>
<!--gt table end-->

<!--/html_preserve-->

    prepped_data %>%
      group_by(`TSS200` = str_detect(UCSC_RefGene_Group, "TSS200") & !is.na(UCSC_RefGene_Group)) %>%
      summarise_at(vars(colon_pwd:endo_pwd), mean, na.rm = TRUE) %>%
      rename(Colon = colon_pwd,
             `Small Intestine` = si_pwd,
             `Endometrium` = endo_pwd) %>%
      mutate(`TSS200` = factor(`TSS200`, c(TRUE, FALSE), c("Yes", "No"))) %>%
      arrange(rev(rownames(.))) %>%
      gt() %>%
      fmt_number(columns = vars(Colon, `Small Intestine`, `Endometrium`), decimals = 3) %>%
      tab_spanner(
        label = "Manhattan distance",
        columns = vars(Colon, `Small Intestine`, `Endometrium`)
      )

<!--html_preserve-->
<style>html {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
}

#svdhjkbuov .gt_table {
  display: table;
  border-collapse: collapse;
  margin-left: auto;
  margin-right: auto;
  color: #000000;
  font-size: 16px;
  background-color: #FFFFFF;
  /* table.background.color */
  width: auto;
  /* table.width */
  border-top-style: solid;
  /* table.border.top.style */
  border-top-width: 2px;
  /* table.border.top.width */
  border-top-color: #A8A8A8;
  /* table.border.top.color */
}

#svdhjkbuov .gt_heading {
  background-color: #FFFFFF;
  /* heading.background.color */
  border-bottom-color: #FFFFFF;
}

#svdhjkbuov .gt_title {
  color: #000000;
  font-size: 125%;
  /* heading.title.font.size */
  padding-top: 4px;
  /* heading.top.padding */
  padding-bottom: 1px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#svdhjkbuov .gt_subtitle {
  color: #000000;
  font-size: 85%;
  /* heading.subtitle.font.size */
  padding-top: 1px;
  padding-bottom: 4px;
  /* heading.bottom.padding */
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#svdhjkbuov .gt_bottom_border {
  border-bottom-style: solid;
  /* heading.border.bottom.style */
  border-bottom-width: 2px;
  /* heading.border.bottom.width */
  border-bottom-color: #A8A8A8;
  /* heading.border.bottom.color */
}

#svdhjkbuov .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  padding-top: 4px;
  padding-bottom: 4px;
}

#svdhjkbuov .gt_col_heading {
  color: #000000;
  background-color: #FFFFFF;
  /* column_labels.background.color */
  font-size: 16px;
  /* column_labels.font.size */
  font-weight: initial;
  /* column_labels.font.weight */
  vertical-align: middle;
  padding: 10px;
  margin: 10px;
}

#svdhjkbuov .gt_sep_right {
  border-right: 5px solid #FFFFFF;
}

#svdhjkbuov .gt_group_heading {
  padding: 8px;
  color: #000000;
  background-color: #FFFFFF;
  /* row_group.background.color */
  font-size: 16px;
  /* row_group.font.size */
  font-weight: initial;
  /* row_group.font.weight */
  border-top-style: solid;
  /* row_group.border.top.style */
  border-top-width: 2px;
  /* row_group.border.top.width */
  border-top-color: #A8A8A8;
  /* row_group.border.top.color */
  border-bottom-style: solid;
  /* row_group.border.bottom.style */
  border-bottom-width: 2px;
  /* row_group.border.bottom.width */
  border-bottom-color: #A8A8A8;
  /* row_group.border.bottom.color */
  vertical-align: middle;
}

#svdhjkbuov .gt_empty_group_heading {
  padding: 0.5px;
  color: #000000;
  background-color: #FFFFFF;
  /* row_group.background.color */
  font-size: 16px;
  /* row_group.font.size */
  font-weight: initial;
  /* row_group.font.weight */
  border-top-style: solid;
  /* row_group.border.top.style */
  border-top-width: 2px;
  /* row_group.border.top.width */
  border-top-color: #A8A8A8;
  /* row_group.border.top.color */
  border-bottom-style: solid;
  /* row_group.border.bottom.style */
  border-bottom-width: 2px;
  /* row_group.border.bottom.width */
  border-bottom-color: #A8A8A8;
  /* row_group.border.bottom.color */
  vertical-align: middle;
}

#svdhjkbuov .gt_striped {
  background-color: #f2f2f2;
}

#svdhjkbuov .gt_from_md > :first-child {
  margin-top: 0;
}

#svdhjkbuov .gt_from_md > :last-child {
  margin-bottom: 0;
}

#svdhjkbuov .gt_row {
  padding: 8px;
  /* row.padding */
  margin: 10px;
  vertical-align: middle;
}

#svdhjkbuov .gt_stub {
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #A8A8A8;
  padding-left: 12px;
}

#svdhjkbuov .gt_summary_row {
  color: #000000;
  background-color: #FFFFFF;
  /* summary_row.background.color */
  padding: 8px;
  /* summary_row.padding */
  text-transform: inherit;
  /* summary_row.text_transform */
}

#svdhjkbuov .gt_grand_summary_row {
  color: #000000;
  background-color: #FFFFFF;
  /* grand_summary_row.background.color */
  padding: 8px;
  /* grand_summary_row.padding */
  text-transform: inherit;
  /* grand_summary_row.text_transform */
}

#svdhjkbuov .gt_first_summary_row {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
}

#svdhjkbuov .gt_first_grand_summary_row {
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #A8A8A8;
}

#svdhjkbuov .gt_table_body {
  border-top-style: solid;
  /* table_body.border.top.style */
  border-top-width: 2px;
  /* table_body.border.top.width */
  border-top-color: #A8A8A8;
  /* table_body.border.top.color */
  border-bottom-style: solid;
  /* table_body.border.bottom.style */
  border-bottom-width: 2px;
  /* table_body.border.bottom.width */
  border-bottom-color: #A8A8A8;
  /* table_body.border.bottom.color */
}

#svdhjkbuov .gt_footnote {
  font-size: 90%;
  /* footnote.font.size */
  padding: 4px;
  /* footnote.padding */
}

#svdhjkbuov .gt_sourcenote {
  font-size: 90%;
  /* sourcenote.font.size */
  padding: 4px;
  /* sourcenote.padding */
}

#svdhjkbuov .gt_center {
  text-align: center;
}

#svdhjkbuov .gt_left {
  text-align: left;
}

#svdhjkbuov .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#svdhjkbuov .gt_font_normal {
  font-weight: normal;
}

#svdhjkbuov .gt_font_bold {
  font-weight: bold;
}

#svdhjkbuov .gt_font_italic {
  font-style: italic;
}

#svdhjkbuov .gt_super {
  font-size: 65%;
}

#svdhjkbuov .gt_footnote_glyph {
  font-style: italic;
  font-size: 65%;
}
</style>

<!--gt table start-->
<table class="gt_table">
<tr>
<th class="gt_col_heading gt_center" rowspan="2" colspan="1">
TSS200
</th>
<th class="gt_col_heading gt_column_spanner gt_center" rowspan="1" colspan="3">
Manhattan distance
</th>
</tr>
<tr>
<th class="gt_col_heading gt_right" rowspan="1" colspan="1">
Colon
</th>
<th class="gt_col_heading gt_right" rowspan="1" colspan="1">
Small Intestine
</th>
<th class="gt_col_heading gt_NA" rowspan="1" colspan="1">
Endometrium
</th>
</tr>
<tbody class="gt_table_body">
<tr>
<td class="gt_row gt_center">
Yes
</td>
<td class="gt_row gt_right">
0.055
</td>
<td class="gt_row gt_right">
0.063
</td>
<td class="gt_row gt_right">
0.051
</td>
</tr>
<tr>
<td class="gt_row gt_center gt_striped">
No
</td>
<td class="gt_row gt_right gt_striped">
0.101
</td>
<td class="gt_row gt_right gt_striped">
0.115
</td>
<td class="gt_row gt_right gt_striped">
0.094
</td>
</tr>
</tbody>
</table>
<!--gt table end-->

<!--/html_preserve-->

Figure 1
--------

    distance_length <- gened_data %>%
      filter(!is.na(UCSC_RefGene_Name)) %>%
      select(CHR, MAPINFO, UCSC_RefGene_Name, colon_pwd:endo_pwd) %>%
      group_by(UCSC_RefGene_Name, CHR) %>%
      mutate(location = MAPINFO - MAPINFO[1]) %>%
      mutate(location_bin = cut(location, c(-1, seq(0, 20000, by = 100)), include.lowest = TRUE)) %>%
      group_by(location_bin) %>%
      summarise(mean_colon = mean(colon_pwd),
                mean_si = mean(si_pwd),
                mean_endo = mean(endo_pwd),
                n_obs = n()) %>%
      mutate(location_bin = as.numeric(str_extract(location_bin, "(?<=,)(.*)(?=])")))

    ## Warning: Factor `location_bin` contains implicit NA, consider using
    ## `forcats::fct_explicit_na`

    distance_length %>%
      drop_na() %>%
      pivot_longer(mean_colon:mean_endo) %>% 
      mutate(name = factor(name, c("mean_colon", "mean_si", "mean_endo"),
                                 c("Colon", "Small Intestine", "Endometrium"))) %>%
      sample_frac() %>%
      ggplot(aes(location_bin, value, color = name, fill = name)) +
      geom_point(key_glyph = draw_key_rect) +
      ggforce::facet_zoom(location_bin < 5000) +
      theme_light() +
      theme(legend.position = "bottom") +
      labs(title = NULL,
           fill = NULL,
           color = NULL,
           x = "Distance by hg19 position from first cpg site (binned to 100s)",
           y = "Average Manhattan Distance", 
           caption = "Figure 1: Average Manhattan distance for single CpGs as a function of position relative to gene.") +
      scale_color_brewer(palette = "Set2") +
      scale_fill_brewer(palette = "Set2")

![](analysis-script_files/figure-markdown_strict/figure1-1.png)

    ggsave(filename = here("figures", "figure-1.png"), width = 7, dpi = 320, height = 4.326)

Figure 2
--------

    permed_data %>%
      pivot_longer(dplyr::matches("1|3")) %>%
      separate(name, into = c("type", "method"), sep = "_v") %>%
      mutate(method = factor(method, c("1", "3"),
                                 c("Naive bootstrap",
                                   "Adjusted bootstrap")),
             type = factor(type, c("colon_pwd", "si_pwd", "endo_pwd"),
                                 c("Colon", "Small Intestine", "Endometrium"))) %>%
      ggplot(aes(value, fill = type)) +
      geom_histogram(bins = 50, color = "grey30") +
      facet_grid(method~type, scales = "free") +
      theme_minimal() +
      labs(x = "p-value",
           y = "Number of genes",
           title = NULL,
           caption = "Figure 2: Distribution of boot-strapping p-values for genes as a function of tissue type and method used.")  +
      theme(legend.position = "none") +
      scale_color_brewer(palette = "Set2") +
      scale_fill_brewer(palette = "Set2")

    ## Warning: Removed 59350 rows containing non-finite values (stat_bin).

![](analysis-script_files/figure-markdown_strict/figure2-1.png)

    ggsave(filename = here("figures", "figure-2.png"), width = 7, dpi = 320, height = 4.326)

    ## Warning: Removed 59350 rows containing non-finite values (stat_bin).

Figure 3
--------

    c(colon = length(enrich_colon$Description),
      si = length(enrich_si$Description),
      endo = length(enrich_endo$Description)) %>%
      sort()

    ## colon    si  endo 
    ##   140   141   158

    list(
      Colon             = enrich_colon$Description,
      `Small Intestine` = enrich_si$Description,
      Endometrium       = enrich_endo$Description
    ) %>%
      UpSetR::fromList() %>%
      UpSetR::upset(sets.bar.color = RColorBrewer::brewer.pal(3, "Set2")[c(3, 2, 1)])
    grid::grid.text("Figure 3: Frequency of overlap between pathways that are called \nas most conserved for each tissue type.",x = 0.65, y=0.05, gp=grid::gpar(fontsize=10))

![](analysis-script_files/figure-markdown_strict/figure3-1.png)

    png(file=here("figures", "figure-3.png"), width = 7, height = 4.326, units = "in", res = 320) # or other device
    list(
      Colon             = enrich_colon$Description,
      `Small Intestine` = enrich_si$Description,
      Endometrium       = enrich_endo$Description
    ) %>%
      UpSetR::fromList() %>%
      UpSetR::upset(sets.bar.color = RColorBrewer::brewer.pal(3, "Set2")[c(3, 2, 1)])
    grid::grid.text("Figure 3: Frequency of overlap between pathways that are called \nas most conserved for each tissue type.",x = 0.65, y=0.05, gp=grid::gpar(fontsize=10))
    dev.off()

    ## quartz_off_screen 
    ##                 2

Figure 4
--------

    ggg <- enrich_si %>%
      ReactomePA::emapplot(showCategory = 30, layout = "graphopt") +
      labs(caption = "Figure 4: Relationship between pathways that are called as conserved in small intestine tissue.",
           title = NULL)

    ggg$layers <- ggg$layers[1:2]

    ggg + 
      ggraph::geom_node_label(aes_(label = ~name), repel = TRUE, size = 2) +
      scale_color_gradientn(colors = "black") +
      guides(color = "none",
             size = "none")

    ## Scale for 'colour' is already present. Adding another scale for
    ## 'colour', which will replace the existing scale.

![](analysis-script_files/figure-markdown_strict/figure4-1.png)

    ggsave(filename = here("figures", "figure-4.png"), width = 7, dpi = 320, height = 4.326)

Figure 5
--------

    figure5_data <- permed_data %>% 
      left_join(by = c("UCSC_RefGene_Name" = "gene"),
      promoted_data
      ) %>%
      left_join(expression_atlas %>% dplyr::select(`Gene Name`, colon, endometrium, `small intestine`),
                by = c("UCSC_RefGene_Name" = "Gene Name")) 

    bind_rows(
      figure5_data %>%
      select(colon_pwd_v3, colon_pro, colon) %>%
      pivot_longer(-colon) %>%
      drop_na() %>%
      rename(expression = colon),
      figure5_data %>%
      select(si_pwd_v3, si_pro, `small intestine`) %>%
      pivot_longer(-`small intestine`) %>%
      drop_na() %>%
      rename(expression = `small intestine`),
      figure5_data %>%
      select(endo_pwd_v3, endo_pro, endometrium) %>%
      pivot_longer(-endometrium) %>%
      drop_na() %>%
      rename(expression = endometrium)
    ) %>%
      group_by(name) %>%
      mutate(value = factor(as.numeric(cut_number(value, 10)))) %>%
      separate(name, c("type", "method"), sep = "_", extra = "merge") %>%
      mutate(type = factor(type, c("colon", "si", "endo"), c("Colon", "Small Intenstine", "Endometrium")),
             method = factor(method, c("pwd_v3", "pro"), c("Boot-strapped values", "Promoter region"))) %>%
      ggplot(aes(value, expression, fill = type)) +
      geom_boxplot() +
      scale_y_log10() +
      guides(color = "none") +
      facet_grid(method ~ type) +
      labs(x = "bin", y = "Expression",
           title = NULL,
           caption = "Figure 5: Relationship between conservation and expression according to tissue type.") +
      theme_minimal() +
      scale_fill_brewer(palette = "Set2") +
      theme(legend.position = "none")

![](analysis-script_files/figure-markdown_strict/figure5-1.png)

    ggsave(filename = here("figures", "figure-5.png"), width = 7, dpi = 320, height = 4.326)
