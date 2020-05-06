import React, { Component } from 'react';
import './App.css';
import Landing from './Landing';
import Footer from './Footer';
import Station from './Station';
import 'mdbreact';


export default class App extends Component {
	constructor(props) {
		super(props);
		this.state = {
			currentDateTime: new Date(),
			program_version: '',
			web_version: '',
		};

		this.get_versions = this.get_versions.bind(this);
		this.get_versions();
	}

	get_versions() {
		fetch('/version')
			.then(res => res.json())
			.then((res) => {
				this.setState({ ...res })
			});
	}

	render() {
		return (
			<div>
				<Landing />
				<Station />
				<Footer { ...this.state } />
			</div>
		);
	}
}
