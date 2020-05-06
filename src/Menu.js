import React, { Component } from 'react';
import { Image } from 'cloudinary-react';
import { Link } from 'react-router-dom';
import { PropTypes } from 'prop-types';
import 'mdbreact';


export default class Menu extends Component {
	static propTypes = {
		stations: PropTypes.arrayOf(
			PropTypes.shape({
				key: PropTypes.number.isRequired,
				name: PropTypes.string.isRequired,
				publicId: PropTypes.string.isRequired,
				href: PropTypes.string.isRequired,
				color: PropTypes.string.isRequired,
			})
		)
	};

	render() {
		return (
			<div className='row row-cols-1 row-cols-md-3'>
				{this.props.stations.map((station, i) =>
					<div className='col' key={i}>
						<div className='card mb-4'>
							<div className='view overlay'>
								<Image className='img-fluid' cloudName='kdphotography-assets' publicId={station.publicId} dpr='auto' responsive width='auto' crop='scale' responsiveUseBreakpoints='true' />
								<Link to={station.href} onClick={(e) => this.props.changeView(e, i+1)}>
									<div className='mask rgba-black-slight'></div>
								</Link>
							</div>
							<div className='card-body text-center'>
								<Link
									id={station.href.slice(1)}
									to={station.href}
									className={`btn btn-${station.color} btn-md`}
									onClick={(e) => this.props.changeView(e, i+1)}
									>
									{station.name}
								</Link>
							</div>
						</div>
					</div>
				)}
			</div>
		);
	}
}
