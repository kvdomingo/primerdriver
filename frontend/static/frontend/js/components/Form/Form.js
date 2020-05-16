import React, { Component } from 'react';
import { Link } from 'react-router-dom';
import {
    MDBContainer as Container,
    MDBTypography as Typography,
    MDBIcon as Icon,
} from 'mdbreact';


export default class Form extends Component {
    render() {
        return (
            <Container className='mb-5'>
                <Link className='btn btn-blue-grey btn-rounded mb-4' to='/' id='back'>
					<Icon fas icon='arrow-left' className='mr-2' /> main menu
				</Link>
				<Typography tag='h2' variant='h2-responsive' className='text-md-center py-2 ml-md-4 d-md-inline'>
					{this.props.title}
				</Typography>
                <form
                    id='form'
                    onChange={this.props.handleValidate}
                    onMouseUp={this.props.handleValidate}
                    onKeyUp={this.props.handleValidate}
                    onSubmit={this.props.handleSubmit}
                    >

                    {this.props.children}

                    <div className='text-center my-3'>
    					<input
    						type='reset'
    						id='reset'
    						onClick={this.props.handleReset}
    						className='btn btn-warning text-dark'
    						value='Reset'
    						/>
                        <input
    						type='submit'
    						id='submit'
                            onClick={this.props.handleSubmit}
    						className='btn btn-primary'
    						value='Submit'
                            disabled={!this.props.isValid}
    						/>
    				</div>
                </form>
            </Container>
        );
    }
}
