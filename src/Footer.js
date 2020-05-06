import React, { Component } from 'react';
import { Image } from 'cloudinary-react';
import 'mdbreact';


export default class Footer extends Component {
	constructor(props) {
		super(props)
		this.state = {
			copyYear: this.props.currentDateTime.getFullYear(),
		}
	}

	render() {
		return (
			<footer className="page-footer font-small bg-dark">
		        <div className="container text-center text-md-left pt-5 pb-4">
		            <div className="row">
		                <div className="col-md-6 mt-md-0 mt-3">
		                    <a href="/">
								<Image className='pb-1' cloudName='kdphotography-assets' publicId='primerdriver/PrimerDriver_logo' height={90} />
		                    </a>
		                    <p>Automated design of mutagenic PCR primers</p>
		                </div>
		                <div className="col-md-3 mb-md-0 mb-3"></div>
		                <hr className="clearfix w-100 d-md-none pb-3" />
		                <div className="col-md-3 mb-md-0 mb-3">
		                    <h5 className="text-uppercase">Links</h5>
		                    <ul className="list-unstyled">
		                        <li>
		                            <a href="https://github.com/kvdomingo/primerdriver/releases" target="_blank" rel="noopener noreferrer">Download</a>
		                        </li>
		                        <li>
		                            <a href="https://kvdomingo.github.io/primerdriver/" target="_blank" rel="noopener noreferrer">Documentation</a>
		                        </li>
		                    </ul>
		                </div>
		            </div>
		        </div>
		        <div className="footer-copyright text-center py-3 px-3">
		            PrimerDriver v{this.props.program_version} {this.props.web_version} &copy; {this.state.copyYear} <a href="mailto:kvdomingo@up.edu.ph">Kenneth Domingo</a> &amp; <a href="mailto:ngutierrez@evc.pshs.edu.ph">Nomer Gutierrez</a>
		        </div>
		    </footer>
		);
	}
}
