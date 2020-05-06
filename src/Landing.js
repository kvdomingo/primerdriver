import React, { Component } from 'react';
import { Image } from 'cloudinary-react';
import 'mdbreact';


export default class Landing extends Component {
	render() {
		return (
			<div className="container text-md-center py-md-5">
				<Image className='img-fluid' cloudName='kdphotography-assets' publicId='primerdriver/PrimerDriver_logo' dpr='auto' responsive width='auto' crop='scale' responsiveUseBreakpoints='true' />
			    <div className="my-4">
			        <p>
			            <b>PrimerDriver</b> ties together key primer design tools and protocols to automate the design of mutagenic PCR primers. This allows input of sequences from different origins (i.e., FASTA, FASTQ, and manual input). The tool can accommodate both DNA & protein sequences to incorporate base pair insertions, deletions, and substitutions as specified by the user.
			        </p>
			        <p>
			            <b>PrimerDriver</b> can design primer pairs and compute for all oligonucleotide sequences that incorporate the desired mutations. This can cater to an array of primer designs from random mutations, site-directed single mutagenesis, batch design site-directed mutagenesis, to multiple-site mutagenesis. Lastly, you can choose to download the results in different file formats - HTML, fasta, and pdf.
			        </p>
			    </div>
			    <a href="#app"><button type="button" className="btn btn-indigo btn-rounded">Get Started</button></a>
			</div>
		);
	}
}
