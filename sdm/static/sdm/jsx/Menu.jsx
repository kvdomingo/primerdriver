class Menu extends React.Component {
    constructor(props) {
        super(props);
        this.pageId = 0;
    }

    render() {
        return (
            <div className='row row-cols-1 row-cols-md-3'>
                {this.props.stations.map((station, i) =>
                    <React.Fragment key={i}>
                        <div className='col'>
                            <div className='card mb-4'>
                                <div className='view overlay'>
                                    <img className='card-img-top' src={station.src} />
                                    <a
                                        href={station.href}
                                        onClick={(e) => this.props.changeView(e, i+1)}
                                    >
                                        <div className='mask rgba-black-slight'></div>
                                    </a>
                                </div>
                                <div className='card-body text-center'>
                                    <a
                                        id={station.href.slice(1)}
                                        href={station.href}
                                        className={`btn btn-${station.color} btn-md`}
                                        onClick={(e) => this.props.changeView(e, i+1)}
                                    >
                                        {station.name}
                                    </a>
                                </div>
                            </div>
                        </div>
                    </React.Fragment>
                )}
            </div>
        );
    }
}
