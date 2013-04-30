package org.caofan.annotate;

import java.util.ArrayList;

import org.apache.commons.lang3.builder.EqualsBuilder;



public final class Reference extends Record {
	//Immutable class
	private final String name;
	private final int cdsStart;
	private final int cdsEnd;
	private final String name2;
	private final int exonCount;
	private final boolean isCoding;
	private final ArrayList<Exon> exons;
	
	public Reference(String name, String chrom, char strand, int start, int end, double score,
			 int cdsStart, int cdsEnd, int exonCount, ArrayList<Exon> exons, String name2, boolean isCoding) {
		super( chrom, start, end, strand, score );
		
		this.name = name;
		this.cdsStart = cdsStart;
		this.cdsEnd = cdsEnd;
		this.exonCount = exonCount;
		this.exons = exons;
		this.name2 = name2;
		this.isCoding = isCoding;
	}
	
	
	public int getCdsStart() {
		return cdsStart;
	}
	
	public int getCdsEnd() {
		return cdsEnd;
	}

	
	public ArrayList<Exon> getExons() {
		return exons;
	}
	
	public String getName() {
		return name;
	}
	
	public String getName2() {
		return name2;
	}
	
	public boolean isCoding() {
		return isCoding;
	}
	
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(name); sb.append("\t");
		sb.append(super.toString());
		sb.append(cdsStart); sb.append("\t");
		sb.append(cdsEnd); sb.append("\t");
		sb.append(exonCount); sb.append("\t");
		for(int i = 0; i < exonCount; i++) {
			sb.append(exons.get(i));
			sb.append(" ");
		}
		sb.append(name2);
		return sb.toString();
		
	}
	
	@Override
	public boolean equals(Object obj) {
		if (obj == null) { return false; }
		if (obj == this) { return true; }
		if (obj.getClass() != getClass()) {
			return false;
		}
		Reference r = (Reference) obj;
		return new EqualsBuilder()
						.appendSuper(super.equals(obj))
						.append(name, r.name)
						.isEquals();
	}
	
}

class Attr {
	final String geneName;
	final String geneType;
	final String txId;
	final String txType;
	
	public Attr(String gn, String gt, String ti, String tt) {
		this.geneName = gn;
		this.geneType = gt;
		this.txId = ti;
		this.txType = tt;
	}
	
	public String getGeneName() {
		return geneName;
	}
	
	public String getGeneType() {
		return geneType;
	}
	
	public String getTxId() {
		return txId;
	}
	
	public String getTxType() {
		return txType;
	}
	
	@Override
	public String toString() {
		return geneName + " " + geneType + " " + txId + " " + txType;
	}
}
