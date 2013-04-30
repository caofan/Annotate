package org.caofan.annotate;

public class QueryRecord extends Record{
	
	final double score;
	final String orig;  //The original line
	
	//These two fields are used to record the annotation information
	String ann;
	String genes;
	
	public QueryRecord(String chrom, int start, int end, char strand, double score, String str) {

		super( chrom, start, end, strand, score );
		this.score = score;
		this.orig = str;
	}
	
	public String getOrig() {
		return orig;
	}
	
	public String getAnn() {
		return ann;
	}
	
	public String getGenes() {
		return genes;
	}
	
	public void setAnn(String ann) {
		this.ann = ann;
	}
	
	public void setGenes(String genes) {
		this.genes = genes;
	}
	
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(orig).append("\t").append(ann).append("\t").append(genes);
		return sb.toString();
	}
}
