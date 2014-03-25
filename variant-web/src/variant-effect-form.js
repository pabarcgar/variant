/*
 * Copyright (c) 2012 Francisco Salavert (ICM-CIPF)
 * Copyright (c) 2012 Ruben Sanchez (ICM-CIPF)
 * Copyright (c) 2012 Ignacio Medina (ICM-CIPF)
 *
 * This file is part of JS Common Libs.
 *
 * JS Common Libs is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * JS Common Libs is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with JS Common Libs. If not, see <http://www.gnu.org/licenses/>.
 */

VariantEffectForm.prototype = new GenericFormPanel();

function VariantEffectForm(args){
    args.analysis = 'hpg-variant.effect';
    GenericFormPanel.prototype.constructor.call(this, args);

	this.id = Utils.genId("VariantEffectForm");
    this.headerWidget =  this.webapp.headerWidget;
    this.opencgaBrowserWidget = this.webapp.headerWidget.opencgaBrowserWidget;
}

VariantEffectForm.prototype.beforeRun = function (){
    //validate regions
    var regions = "";
    var regionPatt = /^([a-zA-Z0-9])+\:([0-9])+\-([0-9])+$/;
    for ( var i = 0; i < this.regionFields.length; i++) {
        var value = this.regionFields[i].getValue();
        if (value!="" && regionPatt.test(value)){
            regions +=value+",";
        }
    }
    if(regions != ''){
        this.paramsWS["region"] = regions.substring(0, regions.length - 1);
    }
    if(this.paramsWS["region"]==''){
        delete this.paramsWS["region"];
    }

    if(this.paramsWS["alleles"]==''){
        delete this.paramsWS["alleles"];
    }
    if(this.paramsWS["coverage"]==''){
        delete this.paramsWS["coverage"];
    }
    this.paramsWS["gene"] = this.paramsWS["gene"].replace(/\s/gm, "").replace(/(^,*)|(,*$)/g, "");
    if(this.paramsWS["gene"]==''){
        delete this.paramsWS["gene"];
    }
    if(this.paramsWS["quality"]==''){
        delete this.paramsWS["quality"];
    }
    if(this.paramsWS["maf"]==''){
        delete this.paramsWS["maf"];
    }

    if(Ext.getCmp("Only SNPs_"+this.id).getValue()){
        this.paramsWS["snp"] = "include";
    }
    if(Ext.getCmp("Only Non-SNPs_"+this.id).getValue()){
        this.paramsWS["snp"] = "exclude";
    }
    if(Ext.getCmp("All_"+this.id).getValue()){
        delete this.paramsWS["snp"];
    }

    /*Input data filter options*/

    /*END Input data filter options*/

    /*Output options*/
    var soTerms = [];
    Ext.getCmp('outputOptions'+this.id).items.each(function(item) {
        if(!item.isDisabled() && item.inputValue != null && item.inputValue != "" && !item.getValue()){
            soTerms.push(item.inputValue);
        }
    });
    if(soTerms.length > 0){
        this.paramsWS["exclude"] = soTerms.toString();
    }

    this.paramsWS["config"] = '/httpd/bioinfo/opencga/analysis/hpg-variant/bin';

    delete this.paramsWS["outputOptions"];
    /*END Output options*/

    if(this.testing){
        console.log("Watch out!!! testing flag is on, so job will not launched.")
    }
};


VariantEffectForm.prototype.getPanels = function (){
	var items = [
	             	this._getSpeciesForm(),
	             	this._getBrowseForm(),
	             	this._getFilterForm(),
	             	this._getOutputForm()
	             ];

	var form = Ext.create('Ext.panel.Panel', {
		margin:"15 0 5 0",
		border:false,
//		layout:{type:'vbox', align: 'stretch'},
		buttonAlign:'center',
		width:"99%",
		//height:900,
		//width: "600",
		items:items,
        defaults:{
            margin:'0 0 15 0'
        }
	});
	
	return [this._getExampleForm(),form];
};
VariantEffectForm.prototype._getSpeciesForm = function (){
	var _this=this;
	
	var checkFlags = function(value){
		var outputOptions = Ext.getCmp('outputOptions'+_this.id);
		if(value!="hsa"){
            outputOptions.child('checkboxfield[inputValue=TF_binding_site_variant]').setValue(false).disable();
            outputOptions.child('checkboxfield[inputValue=miRNA_target_site]').setValue(false).disable();
            outputOptions.child('checkboxfield[id=other_regulatory]').setValue(false).disable();

            outputOptions.child('checkboxfield[inputValue=SNP]').setValue(false).disable();
            outputOptions.child('checkboxfield[id=uniprot_natural_variants]').setValue(false).disable();

            outputOptions.child('checkboxfield[id=phenotypic_annotated_SNPs]').setValue(false).disable();
            outputOptions.child('checkboxfield[id=disease_mutations]').setValue(false).disable();
		}else{
            outputOptions.child('checkboxfield[inputValue=TF_binding_site_variant]').setValue(false).enable();
            outputOptions.child('checkboxfield[inputValue=miRNA_target_site]').setValue(false).enable();
            outputOptions.child('checkboxfield[id=other_regulatory]').setValue(false).enable();

            outputOptions.child('checkboxfield[inputValue=SNP]').setValue(false).enable();
            outputOptions.child('checkboxfield[id=uniprot_natural_variants]').setValue(false).enable();

            outputOptions.child('checkboxfield[id=phenotypic_annotated_SNPs]').setValue(false).enable();
            outputOptions.child('checkboxfield[id=disease_mutations]').setValue(false).enable();
		}
	};

	var speciesForm = Ext.create('Ext.panel.Panel', {
		title:"Species",
        header:this.headerFormConfig,
		border:true,
		padding:"5 0 0 0",
		bodyPadding:10,
		items: []
	});
	
	$.ajax({url:CELLBASE_HOST+"/latest/species?of=json",success:function(data, textStatus, jqXHR){
		// Create the combo box, attached to the states data store
		var objdata = JSON.parse(data);
		for ( var i = 0; i < objdata.length; i++) {
			objdata[i].sciAsembly = objdata[i].scientific+" ("+objdata[i].assembly+")";
		}
		var species = Ext.create('Ext.data.Store', {
			autoLoad: true,
		    fields: ['species', 'common','scientific','assembly','sciAsembly'],
		    data : objdata
		});
		var speciesCombo = Ext.create('Ext.form.field.ComboBox', {
			id:_this.id+"speciesCombo",
		    name:'species',
		    fieldLabel: 'Choose Species',
		    displayField: 'sciAsembly',
		    valueField: 'species',
            editable:false,
		    width:350,
            allowBlank:false,
		    store: species,
			listeners:{
				 change:function(){
					 if(this.getValue()){
						 checkFlags(this.getValue());
			  		}
				 }
			 }
		});
		speciesCombo.select(speciesCombo.getStore().data.items[0]);
		speciesForm.add(speciesCombo);
	},error:function(jqXHR, textStatus, errorThrown){console.log(textStatus);}});
	
  	return speciesForm;
};

VariantEffectForm.prototype._getBrowseForm = function (){
    var _this = this;

    var note1 = Ext.create('Ext.container.Container', {
        html:'<p>Please select a file from your <span class="info">server account</span> using the <span class="emph">Browse</span> button.</p>'
    });


    var formBrowser = Ext.create('Ext.panel.Panel', {
        title:"Select your data",
        header:this.headerFormConfig,
        border:true,
        padding:"5 0 0 0",
        bodyPadding:10,
        items: [note1,this.createOpencgaBrowserCmp({
            fieldLabel: 'Input VCF file:',
            dataParamName: 'vcf-file',
            id: this.id + 'vcf-file',
            mode: 'fileSelection',
            allowedTypes: ['vcf'],
            allowBlank:false
        })]
    });
    return formBrowser;
};


VariantEffectForm.prototype._getExampleForm = function (){
	var _this = this;

	var example1 = Ext.create('Ext.Component', {
		html:'<span class="s140"><span class="btn btn-default">Load</span> &nbsp; VCF file with ~3500 variants example</span>',
		cls:'dedo',
		listeners:{
			afterrender:function(){
				this.getEl().on("click",function(){_this.loadExample1();Ext.example.msg("Example loaded","");});

			}
		}
	});
	var example2 = Ext.create('Ext.Component', {
		html:'<span class="s140"><span class="btn btn-default">Load</span> &nbsp; VCF file with ~5000 variants example</span>',
		cls:'dedo',
		listeners:{
			afterrender:function(){
				this.getEl().on("click",function(){_this.loadExample2();Ext.example.msg("Example loaded","");});
				
			}
		}
	});

    var exampleForm = Ext.create('Ext.container.Container', {
        bodyPadding: 10,
        cls:'bootstrap',
        items: [this.note1, example1],
        defaults: {margin: '5 0 0 0'},
        margin:'0 0 10 0'
    });
	
	return exampleForm;
};


VariantEffectForm.prototype._getFilterForm = function (){
//    By correspondence to an indel (–indel, values are “include/exclude”)
//    By percentage of samples following dominant or recessive inheritance patterns (–inh-dom or –inh-rec, decimal like 0.1)
//    By minimum allele frequency (MAF) (–maf, decimal like 0.01)
//    By percentage of missing values (–missing, decimal like 0.1)
//    By correspondence to a SNP (–snp, values are “include/exclude”)
	var _this=this;
	var items = [];

//    By number of alleles (–alleles)
	var alleles = Ext.create('Ext.form.field.Number', {
		id:this.id+"alleles",
		fieldLabel: 'Alleles',
        name:'alleles',
		width:500,
		minValue:1,
		allowDecimals:false
	});
	items.push(alleles);

//    By coverage (–coverage)
	var coverage = Ext.create('Ext.form.field.Number', {
		id:this.id+"coverage",
		fieldLabel: 'Coverage (min)',
        name:'coverage',
		width:500,
		minValue:0,
		allowDecimals:false
	});
	items.push(coverage);


//    By gene (–gene, comma-separated list)
    var gene = {
        id: this.id + "gene",
        fieldLabel: 'Gene list (csv)',
        xtype: 'textarea',
        name:'gene',
        flex: 1,
        width:500,
        enableKeyEvents: true,
        value: ''
    };
    items.push(gene);

//    By quality (–quality)
	var quality = Ext.create('Ext.form.field.Number', {
		id:this.id+"quality",
		fieldLabel: 'VCF Quality (min)',
        name:'quality',
		width:500,
		minValue:0,
		allowDecimals:false
	});
	items.push(quality);
	var minAlleles = Ext.create('Ext.form.field.Number', {
		id:this.id+"minAlleles",
		fieldLabel: 'Min Alleles Freq (max)',
        name:'maf',
		width:500,
		minValue:0,
		maxValue:1,
		step:0.01,
		decimalPrecision:12,
		allowDecimals:true
	});
	items.push(minAlleles);

	var radioItems = [];
	radioItems.push(this.createRadio("All","snp",true));
	radioItems.push(this.createRadio("Only SNPs","snp"));
	radioItems.push(this.createRadio("Only Non-SNPs","snp"));
	var radioGroup = Ext.create('Ext.form.RadioGroup', {
		fieldLabel: 'SNP',
		width:500,
		items: radioItems
	});
	items.push(radioGroup);


	var formFilterOptions = Ext.create('Ext.form.Panel', {
		title:"Input data filter options",
        header:this.headerFormConfig,
		border:true,
		padding:"5 0 0 0",
		bodyPadding:10,
		items: items
	});

	//regions
//    By region (–region, comma-separated list of regions; –region-file, path to GFF or BED file; –region-type, by feature type of the region, used along –region-file)
	_this.regionFields = [];
	var region = Ext.create('Ext.form.field.Text', {
		id:this.id+"region",
		fieldLabel: 'Region',
        name:'region',
		width:500,
		emptyText:"chr:start-end",
		regex : /^([a-zA-Z0-9])+\:([0-9])+\-([0-9])+$/
	});
	formFilterOptions.insert(4,region);
	_this.regionFields.push(region);

	var button = Ext.create('Ext.button.Button', {
		text:"Add more regions",
		margin:"0 0 15 105",
		handler: function(){
			var reg = Ext.create('Ext.form.field.Text', {
				fieldLabel: 'Region',
                name:'region',
				width:500,
				emptyText:"chr:start-end",
				regex : /^([a-zA-Z0-9])+\:([0-9])+\-([0-9])+$/
			});
			formFilterOptions.insert(4,reg);
			_this.regionFields.push(reg);
		}
	});
	formFilterOptions.insert(5,button);
	
	return formFilterOptions;
};

VariantEffectForm.prototype._getOutputForm = function (){

	var outputOptions = Ext.create('Ext.form.CheckboxGroup', {
		id: 'outputOptions'+this.id,
		columns: 1,
		vertical: true,
        submitValue:false,
		defaults: {margin: '0 0 0 0'},
		items: [
				{ xtype:"label", html:'<span class="emph">Consequence types</span>'},
		        { boxLabel: 'Non-synonymous coding', name: 'outputOptions', inputValue: 'non_synonymous_codon', checked:true },
		        { boxLabel: 'Synonymous coding', name: 'outputOptions', inputValue: 'synonymous_codon' },
		        { boxLabel: 'Splice sites', name: 'outputOptions', inputValue: 'splice_donor_variant,splice_acceptor_variant,splice_region_variant', checked:true },
		        { boxLabel: 'Stop gained/lost', name: 'outputOptions', inputValue: 'stop_gained,stop_lost' },
		        { boxLabel: 'Upstream', name: 'outputOptions', inputValue: '5KB_upstream_variant', checked:true },
		        { boxLabel: 'Downstream', name: 'outputOptions', inputValue: '5KB_downstream_variant', checked:true },
		        { boxLabel: "5' UTR", name: 'outputOptions', inputValue: '5_prime_UTR_variant', checked:true },
		        { boxLabel: "3' UTR", name: 'outputOptions', inputValue: '3_prime_UTR_variant', checked:true },
		        { boxLabel: 'Non-coding RNA', name: 'outputOptions', inputValue: 'pseudogene,nc_transcript_variant,miRNA,lincRNA', checked:true },
		        { boxLabel: 'Intergenic', name: 'outputOptions', inputValue: 'intergenic_variant' },

		        { xtype:"label", html:'<br><span class="emph">Regulatory</span>'},
		        { boxLabel: 'Jaspar TFBS regions', name: 'outputOptions', inputValue: 'TF_binding_site_variant', id: 'TF_binding_site_variant' },
		        { boxLabel: 'miRNA targets', name: 'outputOptions', inputValue: 'miRNA_target_site', id: 'miRNA_target_site' },
		        { boxLabel: 'Other regulatory regions (CTCF, DNaseI, ...)', name: 'outputOptions', inputValue: 'regulatory_region_variant,DNAseI_hypersensitive_site,RNA_polymerase_promoter', id: 'other_regulatory' },


				{ xtype:"label", html:'<br><span class="emph">Variations</span>'},
		        { boxLabel: 'SNPs', name: 'outputOptions', inputValue: 'SNP', id: 'SNP' },
		        { boxLabel: 'Uniprot Natural Variants', name: 'outputOptions', inputValue: '', id: 'uniprot_natural_variants' },//not yet

		        
				{ xtype:"label", html:'<br><span class="emph">Phenotype and diseases</span>'},
		        { boxLabel: 'Phenotypic annotated SNPs', name: 'outputOptions', inputValue: '', id: 'phenotypic_annotated_SNPs' },//not yet --- "no-phenotype"
		        { boxLabel: 'Disease mutations', name: 'outputOptions', inputValue: '', id: 'disease_mutations' }//not yet --- "no-phenotype"
		       ]
	});

	var pan = Ext.create('Ext.form.Panel', {
		title:"Output options",
        header:this.headerFormConfig,
		border:true,
		padding:"5 0 0 0",
		flex:1,
		bodyPadding:10,
		//cls:'panel-border-top',
		items: outputOptions
	});
	return pan;
};


VariantEffectForm.prototype.loadExample1 = function (){

    Ext.getCmp(this.id + 'vcf-file').setText('<span class="emph">Example 1</span>', false);
    Ext.getCmp(this.id + 'vcf-file' + 'hidden').setValue('example_CHB.exon.2010_03.sites.fixed.vcf');


    Ext.getCmp(this.id + 'jobname').setValue("Example vcf 3500");
    Ext.getCmp(this.id + 'jobdescription').setValue("VCF file with ~3500 variants");

    Ext.getCmp("Only SNPs_"+this.id).setValue(true);
//
//
//
//	Ext.getCmp("Only SNPs_"+this.id).setValue(true);
//	this.fileBrowserLabel.setText('<span class="emph">CHB_exon.vcf</span> <span class="info">(server)</span>',false);
//
//	Ext.getCmp("Non-synonymous coding_"+this.id).setValue(true);
//	Ext.getCmp("Synonymous coding_"+this.id).setValue(true);
//	Ext.getCmp("Splice sites_"+this.id).setValue(true);
//	Ext.getCmp("Stop gained/lost_"+this.id).setValue(true);
//	Ext.getCmp("Upstream_"+this.id).setValue(true);
//	Ext.getCmp("Downstream_"+this.id).setValue(true);
//	Ext.getCmp("5' UTR_"+this.id).setValue(true);
//	Ext.getCmp("3' UTR_"+this.id).setValue(false);
//	Ext.getCmp("Non-coding RNA_"+this.id).setValue(true);
//	Ext.getCmp("Intergenic_"+this.id).setValue(false);
//
//	Ext.getCmp("Jaspar TFBS regions_"+this.id).setValue(true);
//	Ext.getCmp("miRNA targets_"+this.id).setValue(true);
//	Ext.getCmp("Other regulatory regions (CTCF, DNaseI, ...)_"+this.id).setValue(false);
//
//	Ext.getCmp("SNPs_"+this.id).setValue(true);
//	Ext.getCmp("Uniprot Natural Variants_"+this.id).setValue(false);
//
//	Ext.getCmp("Phenotypic annotated SNPs_"+this.id).setValue(false);
//	Ext.getCmp("Disease mutations_"+this.id).setValue(false);
//	Ext.getCmp(this.id+"speciesCombo").select(Ext.getCmp(this.id+"speciesCombo").findRecordByValue("hsa"));
//	console.log(this.paramsWS);
////	this.validateRunButton();

};
VariantEffectForm.prototype.loadExample2 = function (){


    Ext.getCmp(this.id + 'vcf-file').setText('<span class="emph">Example 2</span>', false);
    Ext.getCmp(this.id + 'vcf-file' + 'hidden').setValue('example_1000genomes_5000_variants.vcf');


    Ext.getCmp(this.id + 'jobname').setValue("Example vcf 5000");
    Ext.getCmp(this.id + 'jobdescription').setValue("VCF file with ~5000 variants");


//	this.paramsWS["vcf-file-fileid"] = "example2";
//

//
//
//	Ext.getCmp("Only SNPs_"+this.id).setValue(true);
//	this.fileBrowserLabel.setText('<span class="emph">1000genomes_5000_variants.vcf</span> <span class="info">(server)</span>',false);
//
//	Ext.getCmp("Non-synonymous coding_"+this.id).setValue(true);
//	Ext.getCmp("Synonymous coding_"+this.id).setValue(true);
//	Ext.getCmp("Splice sites_"+this.id).setValue(true);
//	Ext.getCmp("Stop gained/lost_"+this.id).setValue(true);
//	Ext.getCmp("Upstream_"+this.id).setValue(true);
//	Ext.getCmp("Downstream_"+this.id).setValue(true);
//	Ext.getCmp("5' UTR_"+this.id).setValue(true);
//	Ext.getCmp("3' UTR_"+this.id).setValue(true);
//	Ext.getCmp("Non-coding RNA_"+this.id).setValue(true);
//	Ext.getCmp("Intergenic_"+this.id).setValue(false);
//
//	Ext.getCmp("Jaspar TFBS regions_"+this.id).setValue(true);
//	Ext.getCmp("miRNA targets_"+this.id).setValue(true);
//	Ext.getCmp("Other regulatory regions (CTCF, DNaseI, ...)_"+this.id).setValue(false);
//
//	Ext.getCmp("SNPs_"+this.id).setValue(true);
//	Ext.getCmp("Uniprot Natural Variants_"+this.id).setValue(false);
//
//	Ext.getCmp("Phenotypic annotated SNPs_"+this.id).setValue(false);
//	Ext.getCmp("Disease mutations_"+this.id).setValue(false);
//	Ext.getCmp(this.id+"speciesCombo").select(Ext.getCmp(this.id+"speciesCombo").findRecordByValue("hsa"));
//	console.log(this.paramsWS);
////	this.validateRunButton();

};


//VariantEffectForm.prototype.validateRunButton = function (){
//	if(this.paramsWS["vcf-file-fileid"] != null && Ext.getCmp("jobNameField_"+this.id).getValue()!=""){
//		this.runButton.enable();
//	}else{
//		this.runButton.disable();
//	}
////	this.runButton.enable();
//};
//VariantEffectForm.prototype.getCheckValue = function (checkbox){
//	if(checkbox.getValue())
//		return null;
//	return "";
//};


//help create functions
//VariantEffectForm.prototype.createTextAreas = function (name, emptyText){
//	var tb = Ext.create('Ext.form.field.TextArea', {
//		id:name+"_"+this.id,
//		fieldLabel : name,
//		name : name,
//		width:500,
//		emptyText:emptyText
////		allowBlank: false
//	});
//	return tb;
//};
//VariantEffectForm.prototype.createTextField = function (name, emptyText){
//	var tb = Ext.create('Ext.form.field.Text', {
//		id:name+"_"+this.id,
//		fieldLabel : name,
//		name : name,
//		width:500,
//		emptyText:emptyText
////		allowBlank: false
//	});
//	return tb;
//};
