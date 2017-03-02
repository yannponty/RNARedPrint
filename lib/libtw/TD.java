import nl.uu.cs.treewidth.*;
import nl.uu.cs.treewidth.ngraph.*;
import nl.uu.cs.treewidth.input.*;
import nl.uu.cs.treewidth.output.*;
import nl.uu.cs.treewidth.algorithm.*;
import nl.uu.cs.treewidth.input.GraphInput.InputData;

import java.io.StringWriter;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.io.IOException;
import java.util.HashMap;

import nl.uu.cs.treewidth.input.GraphInput.InputData;
import nl.uu.cs.treewidth.ngraph.NGraph;
import nl.uu.cs.treewidth.ngraph.NTDBag;
import nl.uu.cs.treewidth.ngraph.NVertex;
import nl.uu.cs.treewidth.ngraph.NEdge;

public class TD{
	
	public static <D> String formatTD( NGraph<NTDBag<D>> g ) {
		
		StringWriter out = new StringWriter();
		
		HashMap<NVertex<NTDBag<D>>,String> bagToName = new HashMap<NVertex<NTDBag<D>>,String>();
		int bagNum = 1;
		for( NVertex<NTDBag<D>> v : g ) {
			String bagName = "bag"+bagNum++;
			bagToName.put( v, bagName );
			out.write( "" + bagName + ": " + v.data.format().replace(" \\n"," ") + "\n" );
		}
		
		out.write( "\n" );
		
		for( NEdge<NTDBag<D>> e : g.edges() ) {
			String nameA = bagToName.get(e.a);
			String nameB = bagToName.get(e.b);
			out.write( "" + nameA + " " + nameB + "\n" );
		}
		
		
		return out.toString();
		
	}
	

	public static void main(String[] argv)
	{
		NGraph<InputData> g = null;

		GraphInput input = new DgfReader(argv[0] );
		try {
			g = input.get();
		} catch( InputException e ) {}
		MaximumMinimumDegreePlusLeastC<InputData> lbAlgo = new MaximumMinimumDegreePlusLeastC<InputData>();
		lbAlgo.setInput( g );
		lbAlgo.run();
		int lowerbound = lbAlgo.getLowerBound();

		GreedyFillIn<InputData> ubAlgo = new GreedyFillIn<InputData>();
		ubAlgo.setInput( g );
		ubAlgo.run();
		int upperbound = ubAlgo.getUpperBound();
		
		NVertexOrder<InputData> permutation = null;

		if( lowerbound == upperbound ) {
			permutation = ubAlgo.getPermutation();
		}
		else {
			QuickBB<InputData> qbbAlgo = new QuickBB<InputData>();
			qbbAlgo.setInput( g );
			qbbAlgo.run();
			permutation = qbbAlgo.getPermutation();
		}
		PermutationToTreeDecomposition<InputData> convertor = new PermutationToTreeDecomposition<InputData>( permutation );
		convertor.setInput( g );
		try{
			convertor.run();
		}
		catch(Exception e){
			System.out.println("Error:");
			e.printStackTrace();
		}
		NGraph<NTDBag<InputData>> decomposition = convertor.getDecomposition();
		try{
			PrintWriter writer = new PrintWriter(argv[1]);
			writer.println(TD.formatTD(decomposition));
			writer.close();
		} catch (IOException e) {
		   // do something
		}		
	}
}
