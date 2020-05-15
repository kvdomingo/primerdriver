import React, { Component } from 'react';
import PropTypes from 'prop-types';


export default class NumberMismatch extends Component {
    static propTypes = {
        handleChange: PropTypes.func.isRequired,
        mismatched_bases: PropTypes.number.isRequired,
    };

    render() {
        return (
            <div className='form-group'>
                <label htmlFor='mismatched_bases'>Number of mismatched bases</label>
                <input
                    className='form-control'
                    min={0}
                    type='number'
                    id='mismatched_bases'
                    name='mismatched_bases'
                    value={this.props.mismatched_bases}
                    onChange={this.props.handleChange}
                    required
                    />
            </div>
        );
    }
}
