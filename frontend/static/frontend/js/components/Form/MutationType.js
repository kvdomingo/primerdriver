import React, { Component } from 'react';
import PropTypes from 'prop-types';


export default class MutationType extends Component {
    static propTypes = {
        handleChange: PropTypes.func.isRequired,
        mutation_type: PropTypes.string.isRequired,
    };

    render() {
        return (
            <div className='form-group'>
                <label htmlFor='mutation_type'>Mutation type</label>
                <select
                    className='browser-default custom-select'
                    id='mutation_type'
                    name='mutation_type'
                    onChange={this.props.handleChange}
                    onMouseUp={(e) => {
                        e.target.name = 'mutation_type';
                        this.props.handleChange(e);
                    }}
                    placeholder='Select mutation type'
                    required
                    value={this.props.mutation_type}
                    >
                    <option name='' value=''>Select mutation type</option>
                    <option value='sub'>Substitution</option>
                    <option value='ins'>Insertion</option>
                    <option value='del'>Deletion</option>
                </select>
            </div>
        );
    }
}
