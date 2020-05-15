import React, { Component } from 'react';


export default class Result extends Component {
    render() {
        return (
            Object.keys(this.props).map((key, i) => (
                <>
                <h2>Henlo</h2>
                <p>{key}: {this.props[key]}</p>
                </>
            ))
        );
    }
}
