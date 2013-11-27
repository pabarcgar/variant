package org.opencb.variant.lib.runners;

import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.opencb.commons.bioformats.pedigree.io.readers.PedDataReader;
import org.opencb.commons.bioformats.variant.VariantStudy;
import org.opencb.commons.bioformats.variant.vcf4.VcfRecord;
import org.opencb.commons.bioformats.variant.vcf4.annotators.VcfGeneNameAnnotator;
import org.opencb.commons.bioformats.variant.vcf4.io.readers.VariantDataReader;
import org.opencb.commons.bioformats.variant.vcf4.io.writers.index.VariantDataWriter;

public class VariantMultiFilterRunner extends VariantRunner {
	
	Map<String, List<VcfRecord>> geneMap;

    public VariantMultiFilterRunner(VariantStudy study, VariantDataReader reader, PedDataReader pedReader, 
            VariantDataWriter writer) {
        this(study, reader, pedReader, writer, null);
    }

    public VariantMultiFilterRunner(VariantStudy study, VariantDataReader reader, PedDataReader pedReader, 
            VariantDataWriter writer,VariantRunner prev) {
        super(study, reader, pedReader, writer, prev);
		this.geneMap = new HashMap<String, List<VcfRecord>>();
    }
	
	@Override
	public List<VcfRecord> apply(List<VcfRecord> batch) throws IOException {
		List<VcfRecord> outputRecords = new LinkedList<VcfRecord>();

		List<String> finishedGenes;
		List<String> geneNames;
		// TODO: en este algoritmo, se comprueban y procesan los finished genes al leer cada variante
		// Podriamos probar si es mas eficiente procesar los finished genes al terminar de leer el batch,
		// aunque consumiría mas ram, si el batch no es excesivamente grande no habria problema
		for (VcfRecord variant : batch) {
			// comprobamos si hay genes que ya han "terminado"
			geneNames = this.variantGeneNames(variant);
			finishedGenes = this.finishedGenes(geneNames, geneMap.keySet());
			if (finishedGenes != null) {
				// filtramos las variantes de los genes terminados y las anadimos a la salida
				outputRecords.addAll(this.applyMultiFilter(finishedGenes, geneMap));				
			}
			
			// if the variant has one of more gene names, add them to the gene map
			if (geneNames != null) {
				this.addVariantAndGenes(variant, geneNames);
			}
		}
		
		return outputRecords;
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
				Arrays.asList(variantInfo.substring(geneNamesTagValuesInitPosition, geneNamesTagEndPosition).split(";"));
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
		// TODO: implementar filtro
		// TODO: no olvidarse de eliminar el gene del map!!
		return passedVariants;
	}

}
