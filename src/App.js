import React, { Component } from 'react';
import './App.css';
import 'mdbreact';


export default class App extends Component {
	constructor(props) {
		super(props)
		this.state = {
			res: null,
			pageId: 0,
			mode: null,
			transitionSpeed: 300,
			transitionName: 'fade',
		}
	}

	changeView = (e, pageId) => {
		e.preventDefault()
		this.setState({ pageId })
	}

	responseCatcher = (e, res, mode) => {
        this.setState({ res: res, mode: mode })
    }

	render() {
		return (
			<div className='container py-3'>
				<a href='!#' className='btn btn-primary'>Primary</a>
			</div>
		);
	}
}
