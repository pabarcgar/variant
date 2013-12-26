package org.opencb.variant.lib.runners;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import org.opencb.commons.bioformats.feature.AllelesCode;
import org.opencb.commons.bioformats.feature.Genotype;
import org.opencb.commons.bioformats.pedigree.Condition;
import org.opencb.commons.bioformats.pedigree.Individual;
import org.opencb.commons.bioformats.pedigree.Pedigree;
import org.opencb.commons.bioformats.pedigree.io.readers.PedDataReader;
import org.opencb.commons.bioformats.variant.VariantStudy;
import org.opencb.commons.bioformats.variant.vcf4.VcfRecord;
import org.opencb.commons.bioformats.variant.vcf4.VcfRecordCoordinateComparator;
import org.opencb.commons.bioformats.variant.vcf4.annotators.VcfGeneNameAnnotator;
import org.opencb.commons.bioformats.variant.vcf4.io.readers.VariantDataReader;
import org.opencb.commons.bioformats.variant.vcf4.io.writers.index.VariantDataWriter;


public class VariantFamiliarGeneFilterRunner extends VariantRunner {
	
	private Map<String, List<VcfRecord>> geneMap;
	
	private Pedigree pedigree;
	
	private Set<String> affectedIndividuals;
	private Set<String> unaffectedIndividuals;
	private Set<String> allIndividuals;

    private Map<VcfRecord,Set<Integer>> sharedAllelesMap;

	private static final int UNAFFECTED = 1;
	
	private static final int AFFECTED = 2;

    public VariantFamiliarGeneFilterRunner(VariantStudy study, VariantDataReader reader, PedDataReader pedReader,
                                           VariantDataWriter writer) {
        this(study, reader, pedReader, writer, null);
    }

    public VariantFamiliarGeneFilterRunner(VariantStudy study, VariantDataReader reader, PedDataReader pedReader,
                                           VariantDataWriter writer, VariantRunner prev) {
        super(study, reader, pedReader, writer, prev);
		this.geneMap = new HashMap<String, List<VcfRecord>>();
		this.pedigree = this.study.getPedigree();
		this.unaffectedIndividuals = this.getIndividualsWithCondition(Condition.UNAFFECTED);
		this.affectedIndividuals = this.getIndividualsWithCondition(Condition.AFFECTED);
		this.allIndividuals = this.pedigree.getIndividuals().keySet();
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
		System.out.println("Total variantes: " + batch.size() + "; Pasan el filtro: " + filteredBatch.size());
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
        System.out.println("Pasan el filtro post: " + lastGenesFilteredRecords.size());
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
		// usamos un set y no una lista para evitar devolver varias veces las mismas variantes
		Set<VcfRecord> passedVariants = new HashSet<VcfRecord>();
		//Map<String, Individual> individuals = this.pedigree.getIndividuals();
		for (String gene : finishedGenes) {
			passedVariants.addAll(this.filterGeneVariants(geneMap.get(gene)));
		}
		// TODO: Aun pueden devolverse variantes repetidas, si estan añadidas en dos genes diferentes, y los dos genes no 
		// estan en la lista de finished. Se devolverian en dos llamadas diferentes a este metodo, a pesar de usar un conjunto aqui
		return passedVariants;
	}
	
	// TODO ¿devolvemos una lista o un conjunto? ¿que es mas eficiente?
	private List<VcfRecord> filterGeneVariants(List<VcfRecord> variantsToBeFiltered) {
		// TODO: cambiar nombre variable y metodo
        // TODO: list o set, y tipo de list, o null
        List<VcfRecord> variantsWhoPassedBothsFilters = new ArrayList<VcfRecord>();

		List<VcfRecord> variantsWhoPassedTheFirstFilter = singleVariantFilter(variantsToBeFiltered);
		if (variantsWhoPassedTheFirstFilter.size() > 1) {
            variantsWhoPassedBothsFilters.addAll(multiVariantFilter(variantsWhoPassedTheFirstFilter));
		}

		return variantsWhoPassedBothsFilters;
	}

    private Set<VcfRecord> multiVariantFilter(List<VcfRecord> variantsToBeFiltered) {
        Set<VcfRecord> res = new HashSet<VcfRecord>();
        List<VcfRecord[]> pairs = variantsPairs(variantsToBeFiltered);
        for (VcfRecord[] pair : pairs) {
            if (!res.contains(pair[0]) || !res.contains(pair[1])) {
                boolean validPair = true;
                for (String unaffectedSample : this.unaffectedIndividuals) {
                    // TODO: cambiar nombre metodo
                    if (sampleHasAllele(unaffectedSample, pair[0], sharedAllelesMap.get(pair[0])) &&
                        sampleHasAllele(unaffectedSample, pair[1], sharedAllelesMap.get(pair[1])))
                    {
                        validPair = false;
                        break;
                    }
                }
                if (validPair) {
                    res.add(pair[0]);
                    res.add(pair[1]);
                }
            }
        }
        // TODO: realmente no es necesario liberar el map, pero, ¿conviene?
        sharedAllelesMap.clear();
        return res;
    }

    private boolean sampleHasAllele(String sample, VcfRecord variant, Set<Integer> alleles) {
        boolean hasAllele = false;
        for (int allele : alleles) {
            Genotype sampleGenotype = variant.getSampleGenotype(sample);
            if (sampleGenotype.getCode() != AllelesCode.ALL_ALLELES_MISSING &&
                (sampleGenotype.getAllele1() == allele || sampleGenotype.getAllele2() == allele))
            {
                // TODO: basta con que el sample tenga alguno de los alelos??
                hasAllele = true;
                break;
            }
        }
        return hasAllele;
    }

    private List<VcfRecord[]> variantsPairs(List<VcfRecord> variantList) {
        // TODO: hacer este metodo mas elegante
        List<VcfRecord[]> pairs = new ArrayList<VcfRecord[]>();
        for (int i=0; i<variantList.size()-1; i ++) {
            for (int j=i+1; j<variantList.size(); j++) {
                VcfRecord[] pair = new VcfRecord[2];
                pair[0] = variantList.get(i);
                pair[1] = variantList.get(j);
                pairs.add(pair);
            }
        }
        return pairs;
    }

	private List<VcfRecord> singleVariantFilter(List<VcfRecord> variantsToBeFiltered) {
		// TODO: ¿array o linked? -> creo que aqui array
		List<VcfRecord> passedVariants = new ArrayList<VcfRecord>();
        this.sharedAllelesMap = new HashMap<VcfRecord, Set<Integer>>();
		//int[] alleles = new int[2];
		for (VcfRecord variant : variantsToBeFiltered) {
            // obtain the shared alternative alleles from the affected samples
			Set<Integer> alternativeAlleles = this.affectedSamplesSharedAlternativeAllelesInVariant(variant);
            if (alternativeAlleles.size() > 0) {
                // check that at least one shared alleles is not homozigous in any unaffected sample
                if (!this.sharedAllelesAreHomozigousInUnaffectedSamples(variant, alternativeAlleles)) {
                    passedVariants.add(variant);
                    this.sharedAllelesMap.put(variant, alternativeAlleles);
                }
            }
		}
		return passedVariants;
	}

    private Set<Integer> affectedSamplesSharedAlternativeAllelesInVariant(VcfRecord variant) {
        Set<Integer> sharedAlternativeAlleles = new HashSet<Integer>();
        boolean firstAffectedSample = true;
        for (String affectedSample : this.affectedIndividuals) {
            Genotype affectedSampleGenotype = variant.getSampleGenotype(affectedSample);
            if (affectedSampleGenotype.getCode() != AllelesCode.ALL_ALLELES_MISSING) {
                //allAffectedAreMissing = false;
                if (firstAffectedSample) {
                    if (affectedSampleGenotype.isAllele1Ref() && affectedSampleGenotype.isAllele2Ref()) {
                        // if an affected sample hasn't at least one alternative allele,
                        // there are no shared alternative alleles
                        sharedAlternativeAlleles.clear();
                        break;
                    }
                    firstAffectedSample = false;
                    // save the alternative alleles
                    if (!affectedSampleGenotype.isAllele1Ref()) {
                        sharedAlternativeAlleles.add(affectedSampleGenotype.getAllele1());
                    }
                    if (!affectedSampleGenotype.isAllele2Ref()) {
                        sharedAlternativeAlleles.add(affectedSampleGenotype.getAllele2());
                    }
                } else {
                    // remove from the alternative alleles set, all the alleles that are not present in the affected sample
                    for (Iterator<Integer> altAllIterator = sharedAlternativeAlleles.iterator(); altAllIterator.hasNext();) {
                        Integer sharedAllele = altAllIterator.next();
                        if (affectedSampleGenotype.getAllele1() != sharedAllele && affectedSampleGenotype.getAllele2() != sharedAllele) {
                            altAllIterator.remove();
                        }
                    }
                    if (sharedAlternativeAlleles.size() == 0) {
                        // there is no common alternative alleles between all the affected samples
                        break;
                    }
                }
            }
        }
        return sharedAlternativeAlleles;
    }

    private boolean sharedAllelesAreHomozigousInUnaffectedSamples(VcfRecord variant, Set<Integer> sharedAlternativeAlleles) {
        boolean allSharedAllelesHomozigousInUnaffected = false;
        // check the unaffected samples
        for (String unaffectedSample : this.unaffectedIndividuals) {
            Genotype unaffectedSampleGenotype = variant.getSampleGenotype(unaffectedSample);
            if ((unaffectedSampleGenotype.getCode() != AllelesCode.ALL_ALLELES_MISSING) &&
                    (!unaffectedSampleGenotype.isAllele1Ref()) &&
                    (unaffectedSampleGenotype.getAllele1() == unaffectedSampleGenotype.getAllele2())) {
                // check if unaffected sample is homozigous for one of the alternative alleles
                for (Iterator<Integer> altAllIterator = sharedAlternativeAlleles.iterator(); altAllIterator.hasNext();) {
                    Integer sharedAllele = altAllIterator.next();
                    if (unaffectedSampleGenotype.getAllele1() == sharedAllele) {
                        altAllIterator.remove();
                    }
                }
                // if the shared alternative alleles set is empty, every shared alleled is homozigous
                // at least in one of the unaffected samples
                if (sharedAlternativeAlleles.size() == 0) {
                    allSharedAllelesHomozigousInUnaffected = true;
                    break;
                }
            }
        }
        return allSharedAllelesHomozigousInUnaffected;
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
