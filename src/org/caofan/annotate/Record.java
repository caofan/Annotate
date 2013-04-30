package org.caofan.annotate;

public class Record implements Comparable< Record >{
	final String chrom;
	final int start;
	final int end;
	final char strand;
	final double score;
	
	public Record( String chrom, int start, int end, char strand, double score ) {
		this.chrom = chrom;
		this.start = start;
		this.end = end;
		this.strand = strand;
		this.score = score;
	}
	
	public String getChrom() {
		return chrom;
	}
	
	public int getStart() {
		return start;
	}
	
	public int getEnd() {
		return end;
	}
	
	public char getStrand() {
		return strand;
	}
	
	public double getScore() {
		return score;
	}
	
	public String toString() {
		StringBuilder sb = new StringBuilder();
		
		sb.append(chrom); sb.append("\t");
		sb.append(strand); sb.append("\t");
		sb.append(start); sb.append("\t");
		sb.append(end); sb.append("\t");
		sb.append(score);
		return sb.toString();
	}
	
	@Override
	public int compareTo(Record r) {
		if (r == null) return -1;
		int score = chrom.compareTo(r.chrom);
		if (score != 0) { return score; }
		else if (start != r.start) { return start - r.start; }
		else if (end != r.end) { return end - r.end; }
		else { return 0; } 
		
	}
}
