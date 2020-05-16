import React, { Component } from 'react';
import './App.css';
import Landing from './Landing';
import Footer from './Footer';
import Station from './Station';


export default class App extends Component {
	constructor(props) {
		super(props);
		this.state = {
			currentDateTime: new Date(),
			program_version: '',
			web_version: '',
		};
	}

	componentDidMount() {
		fetch('/api/version')
			.then(res => res.json())
			.then((res) => {
				this.setState({ ...res })
			});
	}

	render() {
		return (
			<div>
				<Landing />
				<Station id='app' />
				<Footer { ...this.state } />
			</div>
		);
	}
}
