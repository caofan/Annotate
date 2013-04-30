package org.caofan.annotate;

import org.apache.commons.lang3.builder.EqualsBuilder;

public final class Exon {
	final int start;
	final int end;
	final char strand;
	
	public Exon(int start, int end) {
		this(start, end, '+');
	}

	public Exon(int start, int end, char strand) {
		this.start = start;
		this.end = end;
		this.strand = strand;
	}
	
	public int getStart(){
		return start;
	}
	
	public int getEnd() {
		return end;
	}
	
	public char getStrand() {
		return strand;
	}
	
	public String toString() {
		return start + " " + end + " " + strand;
	}
	
	public boolean equals(Object obj) {
		if (obj == null) { return false; }
		if (obj == this) { return true; }
		if (obj.getClass() != getClass()) {
			return false;
		}
		Exon r = (Exon) obj;
		return new EqualsBuilder()
						.appendSuper(super.equals(obj))
						.append(start, r.start)
						.append(end, r.end)
						.append(strand, r.strand)
						.isEquals();
	}
}
