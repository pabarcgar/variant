package org.opencb.variant.lib.core.formats;

import com.fasterxml.jackson.annotation.JsonProperty;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 * User: aleman
 * Date: 9/10/13
 * Time: 8:35 PM
 * To change this template use File | Settings | File Templates.
 */
public class VariantInfo {

    @JsonProperty
    private String chromosome;
    @JsonProperty
    private int position;
    @JsonProperty
    private String ref;
    @JsonProperty
    private String alt;

    @JsonProperty
    private double stats_maf;
    @JsonProperty
    private double stats_mgf;
    @JsonProperty
    private String stats_allele_maf;
    @JsonProperty
    private String stats_genotype_maf;
    @JsonProperty
    private int stats_miss_allele;
    @JsonProperty
    private int stats_miss_gt;
    @JsonProperty
    private int stats_mendel_err;
    @JsonProperty
    private boolean stats_is_indel;
    @JsonProperty
    private double stats_cases_percent_dominant;
    @JsonProperty
    private double stats_controls_percent_dominant;
    @JsonProperty
    private double stats_cases_percent_recessive;
    @JsonProperty
    private double stats_controls_percent_recessive;
    @JsonProperty
    private String stats_id_snp;

    @JsonProperty
    private List<VariantEffect> effect;

    @JsonProperty
    private HashMap<String,String> genotypes;


    public VariantInfo(String chromosome, int position, String ref, String alt) {
        this.chromosome = chromosome;
        this.position = position;
        this.ref = ref;
        this.alt = alt;

        this.effect = new ArrayList<>();

    }

    public VariantInfo(String chromosome, int position, String ref, String alt, VcfVariantStat stats) {
        this(chromosome, position, ref, alt);
        this.addStats(stats);



    }

    public void addStats(VcfVariantStat stat) {

        this.stats_maf = stat.getMaf();
        this.stats_mgf = stat.getMgf();
        this.stats_allele_maf = stat.getMafAllele();
        this.stats_genotype_maf = stat.getMgfAllele();
        this.stats_miss_allele = stat.getMissingAlleles();
        this.stats_miss_gt = stat.getMissingGenotypes();
        this.stats_mendel_err = stat.getMendelinanErrors();
        this.stats_is_indel = stat.getIndel();
        this.stats_cases_percent_dominant = stat.getCasesPercentDominant();
        this.stats_controls_percent_dominant = stat.getControlsPercentDominant();
        this.stats_cases_percent_recessive = stat.getCasesPercentRecessive();
        this.stats_controls_percent_recessive = stat.getControlsPercentRecessive();
        this.stats_id_snp = stat.getId();
    }

    public List<VariantEffect> getEffect() {
        return effect;
    }

    public void setEffect(List<VariantEffect> effect) {
        this.effect = effect;
    }

    public boolean addEffect(VariantEffect ve) {
        return this.effect.add(ve);


    }

    public String getChromosome() {
        return chromosome;
    }

    public void setChromosome(String chromosome) {
        this.chromosome = chromosome;
    }

    public int getPosition() {
        return position;
    }

    public void setPosition(int position) {
        this.position = position;
    }

    public String getRef() {
        return ref;
    }

    public void setRef(String ref) {
        this.ref = ref;
    }

    public String getAlt() {
        return alt;
    }

    public void setAlt(String alt) {
        this.alt = alt;
    }

    public double getStats_maf() {
        return stats_maf;
    }

    public void setStats_maf(double stats_maf) {
        this.stats_maf = stats_maf;
    }

    public double getStats_mgf() {
        return stats_mgf;
    }

    public void setStats_mgf(double stats_mgf) {
        this.stats_mgf = stats_mgf;
    }

    public String getStats_allele_maf() {
        return stats_allele_maf;
    }

    public void setStats_allele_maf(String stats_allele_maf) {
        this.stats_allele_maf = stats_allele_maf;
    }

    public String getStats_genotype_maf() {
        return stats_genotype_maf;
    }

    public void setStats_genotype_maf(String stats_genotype_maf) {
        this.stats_genotype_maf = stats_genotype_maf;
    }

    public int getStats_miss_allele() {
        return stats_miss_allele;
    }

    public void setStats_miss_allele(int stats_miss_allele) {
        this.stats_miss_allele = stats_miss_allele;
    }

    public int getStats_miss_gt() {
        return stats_miss_gt;
    }

    public void setStats_miss_gt(int stats_miss_gt) {
        this.stats_miss_gt = stats_miss_gt;
    }

    public int getStats_mendel_err() {
        return stats_mendel_err;
    }

    public void setStats_mendel_err(int stats_mendel_err) {
        this.stats_mendel_err = stats_mendel_err;
    }

    public boolean isStats_is_indel() {
        return stats_is_indel;
    }

    public void setStats_is_indel(boolean stats_is_indel) {
        this.stats_is_indel = stats_is_indel;
    }

    public double getStats_cases_percent_dominant() {
        return stats_cases_percent_dominant;
    }

    public void setStats_cases_percent_dominant(double stats_cases_percent_dominant) {
        this.stats_cases_percent_dominant = stats_cases_percent_dominant;
    }

    public double getStats_controls_percent_dominant() {
        return stats_controls_percent_dominant;
    }

    public void setStats_controls_percent_dominant(double stats_controls_percent_dominant) {
        this.stats_controls_percent_dominant = stats_controls_percent_dominant;
    }

    public double getStats_cases_percent_recessive() {
        return stats_cases_percent_recessive;
    }

    public void setStats_cases_percent_recessive(double stats_cases_percent_recessive) {
        this.stats_cases_percent_recessive = stats_cases_percent_recessive;
    }

    public double getStats_controls_percent_recessive() {
        return stats_controls_percent_recessive;
    }

    public void setStats_controls_percent_recessive(double stats_controls_percent_recessive) {
        this.stats_controls_percent_recessive = stats_controls_percent_recessive;
    }
}