class Menu extends React.Component {
    constructor(props) {
        super(props)
        if (this.props.currentPage !== 'main') {
            return null;
        }
    }

    goToPage(page) {
        this.setState({ currentPage: page.split('#')[1] });
    }

    render() {
        return (
            <div className='jumbotron blue-grey lighten-5 my-5'>
                <div className='card-deck'>
                    {this.props.stations.map((station, i) =>
                        <React.Fragment>
                            <div className='card mb-4'>
                                <div className='view overlay'>
                                    <img className='card-img-top cld-responsive' data-src={station.src} />
                                    <a
                                        href={station.href}
                                        onClick={this.goToPage.bind(this, station.href)}
                                    >
                                        <div className='mask rgba-black-slight'></div>
                                    </a>
                                </div>
                                <div className='card-body text-center'>
                                    <a
                                        href={station.href}
                                        className={`btn btn-${station.color} btn-md`}
                                        onClick={this.goToPage.bind(this, station.href)}
                                    >
                                        {station.name}
                                    </a>
                                </div>
                            </div>
                        </React.Fragment>
                    )}
                </div>
            </div>
        );
    }
}
