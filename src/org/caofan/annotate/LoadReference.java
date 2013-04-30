package org.caofan.annotate;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.StringTokenizer;

public class LoadReference {
	public static HashMap<String, ArrayList<Record>> loadReference(String filename) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(filename));
			HashMap<String, ArrayList<Record>> all = new HashMap<String, ArrayList<Record>>();
		
			StringBuilder sb = new StringBuilder();
			String line = br.readLine().trim();
			while(line != null && !line.equals("")) {
				StringTokenizer st = new StringTokenizer(line,"\t");
				st.nextToken(); //The bin field.
				String name = st.nextToken();
				String chrom = st.nextToken();
				char strand = st.nextToken().charAt(0);
				int start = Integer.parseInt(st.nextToken());
				int end = Integer.parseInt(st.nextToken());
				int cdsStart = Integer.parseInt(st.nextToken());
				int cdsEnd = Integer.parseInt(st.nextToken());
				boolean isCoding = true;
				if (cdsStart == cdsEnd) {isCoding = false; cdsStart = start; cdsEnd = end;} //Correct for non-coding transcripts.
				int exonCount = Integer.parseInt(st.nextToken());
				StringTokenizer exonStartSt = new StringTokenizer(st.nextToken(), ",");
				StringTokenizer exonEndSt = new StringTokenizer(st.nextToken(), ",");
				ArrayList<Exon> exons = new ArrayList<Exon>();
				for(int i = 0; i < exonCount; i++) {
					int tempStart = Integer.parseInt(exonStartSt.nextToken());
					int tempEnd = Integer.parseInt(exonEndSt.nextToken());
					
					exons.add(new Exon(tempStart, tempEnd, strand));
				}
				double score = Double.parseDouble(st.nextToken()); //The score field.
				String name2 = st.nextToken();
				Reference newRecord = new Reference(name, chrom, strand, start, end, score, 
						cdsStart, cdsEnd, exonCount, exons, name2, isCoding);
				if(!all.containsKey(chrom)) {
					all.put(chrom, new ArrayList<Record>());
				}
				all.get(chrom).add(newRecord);
				line = br.readLine();
			}
			br.close();
			sortRef(all);
			return all;
		} catch(IOException e) {
			System.err.println(e.getMessage());
			System.exit(1);
			return null;
		}
	}
	
	public static HashMap<String, Attr> loadAttrs(String attrFile) {
		try{
			BufferedReader br = new BufferedReader(new FileReader(attrFile));
			String line = br.readLine();
			HashMap<String, Attr> attrsMap = new HashMap<String, Attr>();
			while(line != null) {
				if (line.startsWith("#")) {
					line = br.readLine();
					continue;
				}
				StringTokenizer st = new StringTokenizer(line);
				st.nextToken(); //geneId
				String geneName = st.nextToken().toUpperCase();
				String geneType = st.nextToken().toUpperCase();
				st.nextToken(); //geneStatus
				String txId = st.nextToken().toUpperCase();
				st.nextToken(); //transcriptName
				String txType = st.nextToken().toUpperCase();
				attrsMap.put(txId, new Attr(geneName, geneType, txId, txType));
				line = br.readLine();
			}
			br.close();
			return attrsMap;
		} catch (IOException e) {
			System.err.println(e.getMessage());
			System.exit(1);
			return null;
		}
	}
	
	public static HashMap< String, ArrayList< Record > > loadRejectRegions( String rejectFile ) {
		
		try {
			HashMap<String, ArrayList<Record>> rejectRegions = new HashMap<String, ArrayList<Record>>();

			BufferedReader br = new BufferedReader(new FileReader(rejectFile));
			String line = br.readLine();

			while (line != null) {
				if (line.startsWith("#")) {
					line = br.readLine();
					continue;
				}
				StringTokenizer st = new StringTokenizer(line.trim());
				String chrom = st.nextToken();
				int start = Integer.parseInt(st.nextToken());
				int end = Integer.parseInt(st.nextToken());
				st.nextToken(); //Name
				double score = Double.parseDouble(st.nextToken());
				char strand = st.nextToken().trim().charAt(0);
				
				if (!rejectRegions.containsKey(chrom)) {
					rejectRegions.put(chrom, new ArrayList<Record>());
				}
				rejectRegions.get(chrom).add(
						new Record(chrom, start, end, strand, score));
				line = br.readLine();
			}
			br.close();
			sortRef( rejectRegions );
			return rejectRegions;
		} catch ( IOException e ) {
			System.err.println(e.getMessage());
			System.exit(1);
			return null;
		}
	}
	
	
	
	private static void sortRef(HashMap<String, ArrayList<Record>> ref) {
		Iterator it = ref.entrySet().iterator();
		while (it.hasNext()) {
			Map.Entry pairs = (Map.Entry)it.next();
			//System.out.println(pairs.getKey());
			ArrayList<Record> a = (ArrayList<Record>) pairs.getValue();
			Collections.sort(a, new RefComparator());
		}
	}
	
	
}

class RefComparator implements Comparator<Record> {
	@Override
	public int compare(Record a, Record b) {
		
		return a.compareTo(b);
	}
}
