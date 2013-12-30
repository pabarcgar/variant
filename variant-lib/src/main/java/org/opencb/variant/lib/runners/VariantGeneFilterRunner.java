package org.opencb.variant.lib.runners;

import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import org.opencb.commons.bioformats.pedigree.Condition;
import org.opencb.commons.bioformats.pedigree.Individual;
import org.opencb.commons.bioformats.pedigree.Pedigree;
import org.opencb.commons.bioformats.pedigree.io.readers.PedDataReader;
import org.opencb.commons.bioformats.variant.VariantStudy;
import org.opencb.commons.bioformats.variant.vcf4.VcfRecord;
import org.opencb.commons.bioformats.variant.vcf4.VcfRecordCoordinateComparator;
import org.opencb.commons.bioformats.variant.vcf4.annotators.VcfGeneNameAnnotator;
import org.opencb.commons.bioformats.variant.vcf4.filters.genefilters.VcfFamiliarGeneFilter;
import org.opencb.commons.bioformats.variant.vcf4.filters.genefilters.VcfGeneFilter;
import org.opencb.commons.bioformats.variant.vcf4.io.readers.VariantDataReader;
import org.opencb.commons.bioformats.variant.vcf4.io.writers.index.VariantDataWriter;


public class VariantGeneFilterRunner extends VariantRunner {
	
	private Map<String, List<VcfRecord>> geneMap;

	private Pedigree pedigree;

    private VcfGeneFilter filter;

	private static final int UNAFFECTED = 1;
	
	private static final int AFFECTED = 2;

    public VariantGeneFilterRunner(VariantStudy study, VcfGeneFilter geneFilter, VariantDataReader reader,
                                   PedDataReader pedReader, VariantDataWriter writer) {
        this(study, geneFilter, reader, pedReader, writer, null);
    }

    public VariantGeneFilterRunner(VariantStudy study, VcfGeneFilter geneFilter, VariantDataReader reader,
                                   PedDataReader pedReader, VariantDataWriter writer, VariantRunner prev) {
        super(study, reader, pedReader, writer, prev);
		this.geneMap = new HashMap<String, List<VcfRecord>>();
		this.pedigree = this.study.getPedigree();
        // TODO: ahora mismo el unico filtro existente lo construimos aqui dentro, pero esto debe actualizarse
        //      y recibir el filtro en el constructor
        // this.filter = filter;
        this.filter = new VcfFamiliarGeneFilter(this.getIndividualsWithCondition(Condition.AFFECTED),
                                                this.getIndividualsWithCondition(Condition.UNAFFECTED));
    }
	
	@Override
	public List<VcfRecord> apply(List<VcfRecord> batch) throws IOException {
		List<VcfRecord> filteredBatch = new LinkedList<VcfRecord>();
        // TODO: cambiar nombre
        Set<VcfRecord> preOutputSet = new TreeSet<VcfRecord>(new VcfRecordCoordinateComparator());
        List<VcfRecord> variantsToRemove = null;
		List<String> finishedGenes;
		List<String> geneNames;
		for (VcfRecord variant : batch) {
			// check if there are "finished" genes
			geneNames = this.variantGeneNames(variant);
			finishedGenes = this.finishedGenes(geneNames, geneMap.keySet());
			if (finishedGenes != null) {
				// filter the variants from finished genes, adding those who pass the filters to the output batch
                preOutputSet.addAll(this.applyMultiFilter(finishedGenes, geneMap));
				// remove the finished genes from the gene map
				for (String finishedGene : finishedGenes) {
					geneMap.remove(finishedGene);
				}
                // add to the output batch (in coordinate order) the variants that are not in "not finished" genes
                boolean found = false;
                for (VcfRecord record : preOutputSet) {
                    for (List<VcfRecord> geneVariants: geneMap.values()) {
                        if (geneVariants.contains(record)) {
                            found = true;
                            break;
                        }
                    }
                    if (!found) {
                        filteredBatch.add(record);
                        if (variantsToRemove == null) {
                            variantsToRemove = new LinkedList<VcfRecord>();
                        }
                        variantsToRemove.add(record);
                    } else {
                        break;
                    }
                }
                // remove the processed variants from the set
                if (variantsToRemove != null) {
                    preOutputSet.removeAll(variantsToRemove);
                    variantsToRemove = null;
                }
			}

			// if the variant has one of more gene names, add them to the gene map
			if (geneNames != null) {
				this.addVariantAndGenes(variant, geneNames);
			}
		}
        // write the filtered batch
        batch.clear();
        if (writer != null) {
            ((VariantDataWriter) writer).writeBatch(filteredBatch);
        }

		return filteredBatch;
	}


    public void post() throws IOException {
        // process the last finished genes
        List<VcfRecord> lastGenesFilteredRecords = new LinkedList<VcfRecord>();
        List<String> lastFinishedGenes = new LinkedList<String>();
        lastFinishedGenes.addAll(geneMap.keySet());
        lastGenesFilteredRecords.addAll(this.applyMultiFilter(lastFinishedGenes, geneMap));
        geneMap.clear();
        // write the last genes filtered records
        if (writer != null) {
            ((VariantDataWriter) writer).writeBatch(lastGenesFilteredRecords);
        }
    }
	
	/**
	 * Add the variant and its genes to the geneMap
	 * @param variant - variant to be added
	 * @param geneNames - List with one or more Gene names
	 */
	private void addVariantAndGenes(VcfRecord variant, List<String> geneNames) {
		List<VcfRecord> variantsFromGene;
		for (String geneName : geneNames) {
			variantsFromGene = geneMap.get(geneName);
			if (variantsFromGene == null) {
				// the gene is not in the map, create a new variant list
				variantsFromGene = new LinkedList<VcfRecord> ();
			}
			variantsFromGene.add(variant);
			geneMap.put(geneName, variantsFromGene);
		}
	}
	
	/**
	 * obtain the gene names list from the variant
	 * @param variant
	 * @return
	 */
	private List<String> variantGeneNames(VcfRecord variant) {
		List<String> geneNames = null;
		
		// the Gene names are contained in a tag in the info field of the VCF record
		String variantInfo = variant.getInfo();
		int geneNamesTagPosition = variantInfo.indexOf(VcfGeneNameAnnotator.GENE_NAMES_TAG);
		if (geneNamesTagPosition > -1) {
			int geneNamesTagValuesInitPosition = geneNamesTagPosition + VcfGeneNameAnnotator.GENE_NAMES_TAG.length() + 1;
			// index of the ';' char that marks the end of the geneNames tag, if exists
			int geneNamesTagEndPosition = variantInfo.indexOf(';', geneNamesTagPosition);
			if (geneNamesTagEndPosition == -1) {
				geneNamesTagEndPosition = variantInfo.length();
			}
			// split the gene names into a String list, using the separator ';'
			geneNames = 
				Arrays.asList(variantInfo.substring(geneNamesTagValuesInitPosition, geneNamesTagEndPosition).split(","));
		}
		
		return geneNames;
	}


	/**
	 * Obtain the genes from the set that are not present in the variant gene name list 
	 * @param geneNamesFromVariant - list of gene names from one variant
	 * @param geneSet - gene name set from gene map
	 * @return List of finished genes
	 */
	private List<String> finishedGenes(List<String> geneNamesFromVariant, Set<String> geneSet) {
		List<String> finishedGenes = null;
		// the gene names of the set that are not present in the variant gene name list are the finished ones
		for (String gene : geneSet) {
			if (geneNamesFromVariant == null || !geneNamesFromVariant.contains(gene)) {
				if (finishedGenes == null) {
					finishedGenes = new LinkedList<String>();
				}
				finishedGenes.add(gene);
			}
		}
		return finishedGenes;
	}
	
	/**
	 * 
	 * @param finishedGenes
	 * @param geneMap
	 * @return
	 */
	private Set<VcfRecord> applyMultiFilter(List<String> finishedGenes, Map<String, List<VcfRecord>> geneMap) {
		Set<VcfRecord> passedVariants = new HashSet<VcfRecord>();
		//Map<String, Individual> individuals = this.pedigree.getIndividuals();
		for (String gene : finishedGenes) {
			passedVariants.addAll(this.filter.apply(geneMap.get(gene)));
		}
		return passedVariants;
	}

	private Set<String> getIndividualsWithCondition(Condition condition) {
		Set<String> filteredIndividualNames = new HashSet<String>();
		Map<String, Individual> individuals = this.pedigree.getIndividuals();
		
		for (String individualName : individuals.keySet()) {
			if (individuals.get(individualName).getCondition() == condition) {
				filteredIndividualNames.add(individualName);
			}
		}
		
		return filteredIndividualNames;
	}

}
